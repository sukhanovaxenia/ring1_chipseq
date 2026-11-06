#!/usr/bin/env python3
"""Search ENCODE bigWigs, then explicitly download selected file accessions."""
import argparse, csv, hashlib, json, re, urllib.parse, urllib.request
from pathlib import Path

BASE = 'https://www.encodeproject.org'

def get_json(url):
    req = urllib.request.Request(url, headers={'Accept':'application/json','User-Agent':'RING1-interpretation/1.0'})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.load(r)

def normalized_assembly(value):
    return {'GRCh38':'hg38','hg38':'hg38','GRCh37':'hg19','hg19':'hg19'}.get(value, value)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('action',choices=['search','download'])
    ap.add_argument('--config',default=str(Path(__file__).resolve().parents[1]/'config.json'))
    ap.add_argument('--query',help='Search terms, e.g. fibroblast H3K27me3')
    ap.add_argument('--limit',type=int,default=25)
    ap.add_argument('--accession',help='One file accession ENCFF..., chosen after metadata inspection')
    ap.add_argument('--label',help='Short descriptive track label')
    a=ap.parse_args()
    cp=Path(a.config).resolve();c=json.loads(cp.read_text());root=(cp.parent/c['project_root']).resolve()
    refs=root/c['output_dir']/'references';refs.mkdir(parents=True,exist_ok=True)
    if a.action=='search':
        if not a.query:ap.error('search requires --query')
        if not 1<=a.limit<=100:ap.error('--limit must be 1..100')
        assembly='GRCh38' if c['genome']=='hg38' else 'hg19'
        params={'type':'File','file_format':'bigWig','assembly':assembly,'status':'released','searchTerm':a.query,'limit':a.limit,'format':'json'}
        url=BASE+'/search/?'+urllib.parse.urlencode(params)
        result=get_json(url)
        (refs/'ENCODE.search.json').write_text(json.dumps({'url':url,'response':result},indent=2))
        fields=['accession','assembly','output_type','dataset','biological_replicates','href','file_size','md5sum']
        with open(refs/'ENCODE.search.tsv','w') as h:
            w=csv.DictWriter(h,fieldnames=fields,delimiter='\t');w.writeheader()
            for item in result.get('@graph',[]):
                row={k:item.get(k,'') for k in fields}
                for k,v in row.items():
                    if isinstance(v,(list,dict)):row[k]=json.dumps(v)
                w.writerow(row)
        print('Search total:',result.get('total'),'Saved first',a.limit,'files in',refs/'ENCODE.search.tsv')
        print('Inspect dataset pages for tissue, treatment, target and output type before downloading.')
        return
    if not a.accession or not re.fullmatch(r'ENCFF[0-9A-Z]+',a.accession):ap.error('Provide an ENCODE file accession with --accession')
    item=get_json(BASE+'/files/'+a.accession+'/?format=json')
    if item.get('file_format')!='bigWig' or item.get('status')!='released':raise ValueError('Expected a released bigWig file')
    if normalized_assembly(item.get('assembly'))!=c['genome']:raise ValueError('Assembly mismatch; no automatic coordinate conversion')
    (refs/(a.accession+'.metadata.json')).write_text(json.dumps(item,indent=2))
    path=refs/(a.accession+'.bw')
    if not path.exists():
        tmp=path.with_suffix('.partial')
        url=urllib.parse.urljoin(BASE,item['href'])
        with urllib.request.urlopen(url,timeout=60) as r,open(tmp,'wb') as h:
            while True:
                chunk=r.read(1024*1024)
                if not chunk:break
                h.write(chunk)
        digest=hashlib.md5()
        with open(tmp,'rb') as h:
            for chunk in iter(lambda:h.read(1024*1024),b''):digest.update(chunk)
        if item.get('md5sum') and digest.hexdigest()!=item['md5sum']:raise ValueError('Download checksum mismatch; partial file retained')
        tmp.replace(path)
    manifest=refs/'external_tracks.tsv'
    fields=['label','assembly','path','source_url','output_type','dataset']
    rows=[]
    if manifest.exists():
        with open(manifest) as h:rows=list(csv.DictReader(h,delimiter='\t'))
    rows=[x for x in rows if x['source_url']!=BASE+'/files/'+a.accession+'/']
    rows.append(dict(label=a.label or a.accession,assembly=c['genome'],path=str(path),source_url=BASE+'/files/'+a.accession+'/',output_type=item.get('output_type',''),dataset=str(item.get('dataset',''))))
    with open(manifest,'w') as h:
        w=csv.DictWriter(h,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
    print('Added external track:',path,'; rerun loci to display it on its own scale.')

if __name__=='__main__':main()
