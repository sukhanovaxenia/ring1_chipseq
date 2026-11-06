#!/usr/bin/env python3
"""Orchestrate an auditable single-end RING1 analysis. Run --help first."""
import argparse, csv, gzip, hashlib, json, os, re, shutil, subprocess, sys, urllib.request
from pathlib import Path

BUNDLE = Path(__file__).resolve().parents[1]

def read_bed(path):
    rows = []
    with open(path) as h:
        for n, line in enumerate(h, 1):
            v = line.split()
            if not v or v[0].startswith('#') or v[0] in ('track', 'browser'):
                continue
            if len(v) < 3:
                raise ValueError(f'{path}:{n}: fewer than 3 fields')
            try:
                start, end = int(v[1]), int(v[2])
            except ValueError:
                raise ValueError(f'{path}:{n}: non-integer BED coordinates')
            if start < 0 or end <= start:
                raise ValueError(f'{path}:{n}: invalid interval')
            rows.append((v[0], start, end))
    return rows

def merge_intervals(rows):
    result = []
    for chrom, start, end in sorted(set(rows)):
        if result and result[-1][0] == chrom and start <= result[-1][2]:
            result[-1] = (chrom, result[-1][1], max(end, result[-1][2]))
        else:
            result.append((chrom, start, end))
    return result

def write_bed(path, rows):
    with open(path, 'w') as h:
        for row in rows:
            h.write('\t'.join(map(str, row)) + '\n')

def run(cmd, stdout=None):
    print('+', ' '.join(map(str, cmd)), flush=True)
    subprocess.run(list(map(str, cmd)), check=True, stdout=stdout)

def capture(cmd):
    return subprocess.check_output(list(map(str, cmd)), text=True).strip()

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('stage', choices=['check','fetch-blacklist','prepare','counts','tracks','annotate','matrices','enrich','loci','report'])
    ap.add_argument('--config', default=str(BUNDLE/'config.json'))
    a = ap.parse_args()
    cp = Path(a.config).resolve()
    c = json.loads(cp.read_text())
    root = (cp.parent/c['project_root']).resolve()
    out = root/c['output_dir']
    if a.stage != 'check':
        for folder in ['logs','beds','bams','tables','figures','tracks','matrices','annotation','loci','references']:
            (out/folder).mkdir(parents=True, exist_ok=True)
    resolve = lambda x: root/x
    if c['genome'] not in ('hg38','hg19'):
        raise ValueError('This bundle supports human UCSC hg38 or hg19 only.')
    if len(c['samples']) != 4 or [s['condition'] for s in c['samples']] != ['BF','BF','HS','HS']:
        raise ValueError('Expected four samples in BF1, BF2, HS1, HS2 order.')
    for s in c['samples']:
        if not re.fullmatch(r'[A-Za-z0-9_]+', s['sample']):
            raise ValueError('Sample labels must contain letters, digits or underscore only.')
    if len({s['sample'] for s in c['samples']}) != 4:
        raise ValueError('Sample labels must be unique.')
    if c['fragment_length'] is not None and int(c['fragment_length']) <= 0:
        raise ValueError('fragment_length must be null or positive.')
    blacklist = resolve(c['blacklist'])
    signature_config = {k:v for k,v in c.items() if k not in ('threads','biological_replicates_confirmed')}
    signature = hashlib.sha256(json.dumps(signature_config,sort_keys=True).encode()).hexdigest()
    state = out/'logs/prepare_config.sha256'
    if a.stage not in ('check','fetch-blacklist','report') and state.exists() and state.read_text().strip()!=signature:
        raise RuntimeError('Analysis configuration changed. Use a NEW output_dir to avoid mixing results.')
    if a.stage not in ('check','fetch-blacklist','prepare','report') and not state.exists():
        raise RuntimeError('Run prepare successfully first.')
    if a.stage == 'fetch-blacklist':
        if blacklist.exists():
            print('Already exists:', blacklist); return
        url = f'https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/{c["genome"]}-blacklist.v2.bed.gz'
        content = gzip.decompress(urllib.request.urlopen(url, timeout=60).read())
        blacklist.parent.mkdir(parents=True, exist_ok=True)
        temp = blacklist.with_suffix('.download'); temp.write_bytes(content)
        read_bed(temp); temp.replace(blacklist)
        (out/'references/blacklist_source.json').write_text(json.dumps({'url':url,'sha256':hashlib.sha256(content).hexdigest()}, indent=2))
        return
    if a.stage == 'check':
        for tool in ['samtools','bedtools','Rscript','bamCoverage','computeMatrix','plotHeatmap','plotProfile']:
            print(tool, shutil.which(tool) or 'MISSING')
        for s in c['samples']:
            for field in ['bam','input','peaks']:
                p=resolve(s[field]); print(s['sample'],field,'OK' if p.is_file() else 'MISSING',p)
        for label,path in {**c['pooled_peaks'],**c['bdgdiff_regions']}.items():
            p=resolve(path); print(label,'OK' if p.is_file() else 'MISSING',p)
        print('blacklist:',blacklist.is_file(),'; biological repeats confirmed:',c['biological_replicates_confirmed'])
        run(['Rscript',BUNDLE/'scripts/install.R','check',c['genome']]); return
    if a.stage in ['counts','annotate','enrich','loci']:
        run(['Rscript',BUNDLE/'scripts/analysis.R',cp,a.stage]); return
    if a.stage == 'prepare':
        if not blacklist.is_file():
            raise FileNotFoundError('Run fetch-blacklist first, or provide the matching blacklist.')
        if state.exists() and state.read_text().strip()!=signature:
            raise RuntimeError('Config changed. Select a NEW output_dir to avoid mixing analysis versions.')
        expected_chr1={'hg38':248956422,'hg19':249250621}[c['genome']]
        allowed=None; sizes=None
        libraries=[]
        for s in c['samples']:
            for kind, field in [('ChIP','bam'),('Input','input')]:
                bam=resolve(s[field]); run(['samtools','quickcheck','-v',bam])
                hdr=capture(['samtools','view','-H',bam])
                d={}
                for line in hdr.splitlines():
                    if line.startswith('@SQ\t'):
                        tags=dict(x.split(':',1) for x in line.split('\t')[1:] if ':' in x)
                        d[tags['SN']]=int(tags['LN'])
                if d.get('chr1')!=expected_chr1:
                    raise ValueError(f'{bam}: chr1 does not match {c["genome"]}; check build and naming.')
                canonical={k:v for k,v in d.items() if re.fullmatch(r'chr([1-9]|1[0-9]|2[0-2]|X|Y)',k)}
                if sizes is None:
                    sizes=canonical; allowed=set(sizes)
                elif canonical!=sizes:
                    raise ValueError('BAM chromosome dictionaries disagree.')
                label=s['sample'] if kind=='ChIP' else s['sample']+'_Input'
                target=out/'bams'/f'{label}.analysis.bam'
                # Explicitly sort after chromosome selection and blacklist filtering.
                if not Path(str(bam)+'.bai').exists() and not bam.with_suffix('.bai').exists() and not Path(str(bam)+'.csi').exists():
                    raise FileNotFoundError(f'Index missing: samtools index {bam}')
                if not target.exists():
                    tmp=out/'bams'/f'{label}.preblacklist.bam'
                    run(['samtools','view','-@',c['threads'],'-b','-q',c['mapq'],'-F','2308','-o',tmp,bam,*sorted(allowed)])
                    part=Path(str(target)+'.partial')
                    with open(part,'wb') as h:
                        run(['bedtools','intersect','-v','-abam',tmp,'-b',blacklist],stdout=h)
                    sorted_part=Path(str(target)+'.sorted.partial')
                    run(['samtools','sort','-@',c['threads'],'-o',sorted_part,part])
                    run(['samtools','quickcheck',sorted_part]); sorted_part.replace(target)
                    part.unlink(); tmp.unlink()
                run(['samtools','index','-@',c['threads'],target])
                libraries.append({'sample':label,'condition':s['condition'],'kind':kind,'source_bam':str(bam),'bam':str(target),'eligible_reads':int(capture(['samtools','view','-c',target]))})
        with open(out/'tables/libraries.tsv','w') as h:
            w=csv.DictWriter(h,fieldnames=list(libraries[0]),delimiter='\t');w.writeheader();w.writerows(libraries)
        with open(out/'references/chrom.sizes','w') as h:
            for k,v in sorted(sizes.items()): h.write(f'{k}\t{v}\n')
        groups=dict(c['pooled_peaks']);groups.update(c['bdgdiff_regions'])
        groups.update({s['sample']:s['peaks'] for s in c['samples']})
        region_rows={};manifest=[]
        for label,path in groups.items():
            p=resolve(path)
            if not p.exists():
                if label in c['pooled_peaks']: raise FileNotFoundError(p)
                print('Optional region set missing:',p); continue
            raw=read_bed(p)
            valid=[r for r in raw if r[0] in allowed]
            if any(r[2]>sizes[r[0]] for r in valid): raise ValueError(f'{p}: coordinates exceed chromosome.')
            temp=out/'beds'/f'{label}.tmp.bed';write_bed(temp,sorted(set(valid)))
            dest=out/'beds'/f'{label}.bed'
            with open(dest,'w') as h:run(['bedtools','intersect','-v','-a',temp,'-b',blacklist],stdout=h)
            temp.unlink(); region_rows[label]=read_bed(dest)
            manifest.append({'group':label,'source':str(p),'bed':str(dest),'original_n':len(raw),'retained_n':len(region_rows[label])})
        union=merge_intervals(region_rows['BF_pool']+region_rows['HS_pool'])
        if not union:raise ValueError('No pooled regions remain after filtering.')
        write_bed(out/'beds/consensus.bed',union)
        with open(out/'tables/region_sets.tsv','w') as h:
            w=csv.DictWriter(h,fieldnames=list(manifest[0]),delimiter='\t');w.writeheader();w.writerows(manifest)
        # Each alignment is counted once with -L, even if it overlaps several intervals.
        for row in libraries:
            in_peaks=int(capture(['samtools','view','-c','-L',out/'beds/consensus.bed',row['bam']]))
            row['reads_in_consensus']=in_peaks
            row['FRiP_consensus_read_overlap']=in_peaks/row['eligible_reads'] if row['eligible_reads'] else 0
        with open(out/'tables/read_QC.tsv','w') as h:
            w=csv.DictWriter(h,fieldnames=list(libraries[0]),delimiter='\t');w.writeheader();w.writerows(libraries)
        state.write_text(signature+'\n'); shutil.copy2(cp,out/'logs/config.used.json'); return
    if a.stage == 'tracks':
        with open(out/'tables/normFactors.txt') as h: rows=list(csv.DictReader(h,delimiter='\t'))
        for row in rows:
            command=['bamCoverage','-b',row['bam'],'-o',out/'tracks'/f'{row["sample"]}.bw','--binSize','10','--normalizeUsing','None','--scaleFactor',row['bamCoverage_scale'],'-p',c['threads']]
            if c['fragment_length'] is not None:command+=['--extendReads',c['fragment_length']]
            run(command)
        return
    if a.stage == 'matrices':
        bws=[out/'tracks'/f'{s["sample"]}.bw' for s in c['samples']]
        if not all(p.exists() for p in bws):raise FileNotFoundError('Run counts and tracks first.')
        jobs=[('consensus',[out/'beds/consensus.bed'],['consensus'],'center')]
        diff=[(k,out/'beds'/f'{k}.bed') for k in c['bdgdiff_regions']]
        diff=[(k,p) for k,p in diff if p.exists() and p.stat().st_size]
        if diff:jobs.append(('bdgdiff',[p for k,p in diff],[k for k,p in diff],'center'))
        if (out/'beds/genes.bed').exists():
            jobs.append(('TSS',[out/'beds/genes.bed'],['genes'],'TSS'))
            jobs.append(('gene_body',[out/'beds/genes.bed'],['genes'],'scaled'))
        sig=[(k,out/'beds'/f'{k}.bed') for k in ['BF_higher','HS_higher']]
        sig=[(k,p) for k,p in sig if p.exists() and p.stat().st_size]
        if sig:jobs.append(('differential',[p for k,p in sig],[k for k,p in sig],'center'))
        for name,beds,labels,ref in jobs:
            matrix=out/'matrices'/f'{name}.gz'
            mode=['scale-regions','--regionBodyLength','5000'] if ref=='scaled' else ['reference-point','--referencePoint',ref]
            run(['computeMatrix',*mode,'-b',c['profile_flank'],'-a',c['profile_flank'],'--binSize',c['profile_bin'],'-R',*beds,'-S',*bws,'--samplesLabel',*[s['sample'] for s in c['samples']],'--sortRegions','keep','-p',c['threads'],'-o',matrix])
            # Shared 98th percentile across ALL tracks; record truncation explicitly.
            import numpy as np
            arrays=[]
            with gzip.open(matrix,'rt') as h:
                for line in h:
                    if line.startswith('@'):continue
                    arrays.append(np.fromstring(line.rstrip().split('\t',6)[6],sep='\t'))
            if not arrays:raise ValueError(f'Empty matrix: {name}')
            vals=np.concatenate(arrays); vals=vals[np.isfinite(vals)]
            if not len(vals):raise ValueError('All matrix values missing; check assembly/chr naming.')
            upper=max(float(np.quantile(vals,0.98)),1e-9)
            (out/'matrices'/f'{name}.scale.json').write_text(json.dumps({'shared_zmax':upper,'rule':'98th percentile of all finite bins across all samples; higher values clipped','regions':labels},indent=2))
            run(['plotHeatmap','-m',matrix,'-out',out/'figures'/f'{name}.heatmap.pdf','--regionsLabel',*labels,'--refPointLabel',ref,'--sortRegions','descend','--sortUsing','mean','--zMin','0','--zMax',upper,'--colorMap','YlOrRd','--missingDataColor','grey','--outFileSortedRegions',out/'matrices'/f'{name}.sorted.bed'])
            run(['plotProfile','-m',matrix,'-out',out/'figures'/f'{name}.profile.pdf','--perGroup','--regionsLabel',*labels,'--refPointLabel',ref,'--plotType','lines','--outFileNameData',out/'matrices'/f'{name}.profile.tsv'])
        return
    if a.stage=='report':
        import html
        groups=[]
        for folder in ['figures','tables','annotation','loci','matrices','logs']:
            links='\n'.join(f'<li><a href="{html.escape(str(p.relative_to(out)),quote=True)}">{html.escape(p.name)}</a></li>' for p in sorted((out/folder).glob('*')) if p.is_file())
            groups.append(f'<h2>{folder}</h2><ul>{links}</ul>')
        (out/'index.html').write_text('<!doctype html><meta charset="utf-8"><title>RING1 results</title><style>body{max-width:1000px;margin:40px auto;font:17px system-ui;line-height:1.6}</style><h1>RING1 BF vs HS</h1><p>logFC &gt; 0: HS выше BF. Результаты требуют проверки качества повторов и происхождения контролей. См. README набора скриптов.</p>'+''.join(groups))
        print(out/'index.html')

if __name__=='__main__':
    main()
