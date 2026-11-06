# RING1 ChIP-seq: фибробласты молочной железы (BF) и крайней плоти (HS)

Анализ связывания RING1 (компонент PRC1) в двух линиях фибробластов человека по данным [GSE40740](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40740) (Pemberton et al., *Genome Biology* 2014, 15:R23):

- **BF** — фибробласты молочной железы взрослой женщины;
- **HS** — неонатальные фибробласты крайней плоти Hs68.

Для каждой линии есть два повтора ChIP RING1 и по контролю (input) к каждому повтору.

**В данном репозитории описывают этапы пострпроцессинга, выполненные в рамках самостоятельной работы** Выравнивание, удаление дубликатов, MACS2 callpeak и MACS2 bdgdiff уже выполнены заранее на сервере в качестве тоже самостоятельного анализа. В репозитории и включены:

- контроль качества повторов;
- общий набор регионов и дифференциальное связывание BF и HS (csaw + edgeR);
- аннотация пиков;
- профили сигнала;
- bdgdiff по отдельным повторам с аннотацией HOMER;
- снимки IGV ключевых локусов.

Описание результатов — в [RESULTS.md](RESULTS.md).

Автор и исполнитель: **Xenia Sukhanova**.

## Структура

```text
.
├── README.md
├── RESULTS.md                основные результаты
├── config.example.json       шаблон конфигурации (скопировать в config.json)
├── environment.yml           окружение conda
├── scripts/
│   ├── run.py                основной конвейер, запуск по этапам
│   ├── analysis.R            csaw + edgeR (дифференциальный анализ), ChIPseeker, GO, Gviz; вызывается из run.py
│   ├── install.R             установка и проверка R-пакетов
│   ├── bdgdiff_reps.sh       MACS2 bdgdiff по отдельным повторам
│   ├── biomart_mapping.R     необязательно: Entrez → Ensembl через BioMart
│   └── encode_refs.py        необязательно: поиск и загрузка треков ENCODE
├── preprocessing/            скрипты предварительных шагов (тримминг, выравнивание, коллинг, bdgdiff), в исходном виде
├── tests/test_helpers.py     тесты вспомогательных функций run.py
└── results/
    ├── tables/               таблицы QC, нормализации и дифференциального анализа
    ├── figures/              графики конвейера (PDF)
    ├── beds/                 наборы регионов (BED3)
    ├── homer/                аннотация bdgdiff по повторам (annotatePeaks.pl)
    ├── tss_profiles/         профили deepTools вокруг TSS для регионов bdgdiff
    ├── igv/                  снимки IGV (PDF)
    ├── png/                  PNG-версии графиков и снимков IGV для RESULTS.md
    └── qc/                   отчёты MultiQC (FastQC до и после тримминга)
```

BAM, bigWig, bedGraph и матрицы deepTools в репозиторий не входят из-за размера.

## Подготовка

Скрипты ожидают каталог проекта (`ring/`) с результатами предыдущих шагов. Репозиторий клонируется внутрь него:

```text
ring/
├── align/chip/dedup/           *.dedup.bam (+ .bai) — 4 ChIP и 4 input
├── peaks/
│   ├── RING1_BF_peaks.narrowPeak, RING1_HS_peaks.narrowPeak   пулированный коллинг
│   ├── diff/                   bdgdiff по пулированным данным
│   └── reps/                   коллинг по повторам (MACS2 -B)
└── ring1-chipseq/              этот репозиторий
```

Предварительные шаги, отсутствующие здесь и предшествующие постпроцессингу:

1. FastQC и тримминг прочтений.
2. Выравнивание на hg38 и удаление дубликатов.
3. MACS2 callpeak:
   - по пулированным повторам каждой линии;
   - по каждому повтору отдельно с ключом `-B`, чтобы получить `*_treat_pileup.bdg` и `*_control_lambda.bdg`.
4. MACS2 bdgdiff BF против HS по пулированным данным.

## Установка

```bash
conda env create -f environment.yml
conda activate ring-interpret
Rscript scripts/install.R install hg38   # доустановит недостающие R-пакеты
Rscript scripts/install.R check hg38
```

Для отдельных шагов нужны также MACS2 (для `bdgdiff_reps.sh`) и HOMER с геномом hg38 (для аннотации).

## Запуск

### 1. Конфигурация

```bash
cp config.example.json config.json
```

Что проверить в `config.json`:

- **Пути.** Все пути задаются относительно `project_root`. По умолчанию это `..`, то есть каталог `ring/`.
- **samples.** Четыре образца строго в порядке BF_rep1, BF_rep2, HS_rep1, HS_rep2. Для каждого указываются ChIP BAM, BAM контроля и narrowPeak повтора.
- **blacklist.** Путь к blacklist ENCODE. Если файла нет, его скачает этап `fetch-blacklist`.
- **Параметры статистики.**
  - `fdr` — порог FDR, по умолчанию 0.05.
  - `minimum_log2_fold_change` — порог TREAT, по умолчанию 1, то есть различие больше чем в 2 раза.
  - `background_bin` — размер окна для нормализации по фону, по умолчанию 10 kb.

Результаты сохраняются в `ring/<output_dir>`. Если после `prepare` буду изменяться параметры, то нужно будут создать (и объявить в конфиге)  новый `output_dir`, так как в ходе каждого запуска происходит проверка и соотнесение параметров конфигурационного файла с локальной структурой.

### 2. Основной пайплнай

Команды запускаются из корня репозитория:

```bash
python scripts/run.py check             # инструменты, пути, R-пакеты
python scripts/run.py fetch-blacklist   # blacklist ENCODE v2, если его ещё нет
python scripts/run.py prepare           # фильтрация BAM (chr1–22, X, Y; MAPQ ≥ 20; без blacklist), наборы регионов BED3, FRiP
python scripts/run.py counts            # csaw: подсчёт в общем наборе регионов, нормализация по фону, edgeR QL + TREAT
python scripts/run.py tracks            # bigWig, нормированные на размер библиотеки и TMM
python scripts/run.py annotate          # ChIPseeker: категории, расстояние до TSS, ближайшие гены
python scripts/run.py matrices          # deepTools: тепловые карты и профили
python scripts/run.py enrich            # GO BP для генов дифференциальных регионов
python scripts/run.py loci              # Gviz: отдельные локусы
python scripts/run.py report            # index.html со ссылками на все результаты
```

У этапа `loci` два режима:

- по умолчанию он берёт до трёх самых значимых регионов в каждом направлении;
- свой список можно задать файлом `ring/interpretation/loci/regions_to_plot.tsv` с колонками `chr`, `start0`, `end` и `label`, разделёнными табуляцией.

### 3. Дифференциальный анализ (csaw + edgeR)

Код находится в `scripts/analysis.R`, в блоке этапа `counts`; отдельно его запускать не нужно. Статистическая часть выполняется, только если в конфиге стоит `"biological_replicates_confirmed": true`. Иначе считаются лишь подсчёты и QC.

1. **Регионы.** Тестируется общий набор `consensus.bed`, то есть объединение пулированных пиков BF и HS после фильтрации по blacklist. Все четыре ChIP-библиотеки оцениваются на одних и тех же координатах.
2. **Подсчёт прочтений.** `csaw::regionCounts` считает прочтения в каждом регионе по отфильтрованным BAM этапа `prepare`. Фрагменты не удлиняются, пока в конфиге `fragment_length: null`.
3. **Нормализация.**
   - Факторы нормализации считаются не по пикам, а по фоновым окнам: `windowCounts` с окном `background_bin` (10 kb), затем `normFactors`. Так различие в силе связывания между линиями не гасится нормализацией.
   - Инпут в модель не входит. Отношение ChIP/input сохраняется только для проверки в `ChIP_Input.exploratory_log2ratio.tsv`.
4. **Фильтрация.** `filterByExpr` оставляет регионы с достаточным числом прочтений в группе (2777 из 4209). Список протестированных регионов записывается в `tested_regions.tsv`.
5. **Модель.** `design = ~ group`, где BF — референсный уровень. `estimateDisp` и `glmQLFit` запускаются с `robust = TRUE`, коэффициент 2 соответствует HS относительно BF, так что `logFC > 0` означает сигнал выше в HS. Повторы с одинаковым номером не считаются парными.
6. **Тест.** Используются два теста:
   - `glmQLFTest` (колонки `PValue` и `FDR`) отвечает на вопрос, отличается ли logFC от нуля;
   - `glmTreat` с `lfc = minimum_log2_fold_change` (колонка `FDR_TREAT`) — превышает ли |log2FC| порог.

   Значимыми считаются регионы с `FDR_TREAT < fdr`. Они делятся на `HS_higher` и `BF_higher` по знаку logFC.

Что получается на выходе:

| Файл | Содержимое |
|---|---|
| `tables/differential.HS_vs_BF.tsv` | все протестированные регионы: logFC, logCPM, F, PValue, FDR, FDR_TREAT, direction |
| `tables/differential.annotated.tsv` | то же с аннотацией ChIPseeker (этап `annotate`) |
| `beds/HS_higher.bed`, `beds/BF_higher.bed` | значимые регионы по направлениям |
| `tables/normFactors.txt`, `tables/normalization_QC.tsv` | размеры библиотек и факторы нормализации; из них же считается масштаб bigWig |
| `figures/dispersion.pdf`, `QL_dispersion.pdf` | BCV и QL-дисперсии |
| `figures/differential.volcano.pdf`, `differential.MA.pdf` | volcano и MA |

### 4. bdgdiff по отдельным повторам

```bash
bash scripts/bdgdiff_reps.sh /path/to/ring
```

Скрипт сравнивает BF_repN и HS_repN, используя pileup и lambda из `peaks/reps/`. Глубину каждой пары он берёт из `*_peaks.xls`, которая соответствует наименьшему числу прочтений ChIP и контроля после фильтрации. На выходе получаются файлы `peaks/reps/diff/RING1_repN_{BF_unique,HS_unique,BF_HS_common}.bed`.

### 5. Аннотация HOMER

```bash
mkdir -p peaks/reps/stats
for f in peaks/reps/diff/*.bed; do
    b=$(basename "$f" .bed)
    annotatePeaks.pl "$f" hg38 -cpu 6 -annStats peaks/reps/stats/$b.stats.txt > annotation/reps/$b.annotated.bed
done
```

### 6. Профили вокруг TSS (deepTools)

Профили в `results/tss_profiles/` построены скриптом `preprocessing/visualize_chip.sh` по регионам пулированного bdgdiff, ±1 kb от TSS:

```bash
computeMatrix reference-point --referencePoint TSS -b 1000 -a 1000 --skipZeros \
    -R peaks/diff/bed3/RING1_BF_unique.bed \
    -S normalization/bin/RING1_BF_rep*.binned.TMM.bigwig -o BF_TSS.gz
plotProfile -m BF_TSS.gz -out TSS_BF_profile.png --perGroup
plotHeatmap -m BF_TSS.gz -out TSS_BF_profile-heatmap.png
```

Для HS команды те же.

### 7. Необязательные шаги

```bash
# Ensembl ID и описания генов; второй аргумент — фиксированный релиз Ensembl
Rscript scripts/biomart_mapping.R config.json 113

# Внешние треки ENCODE: сначала поиск, затем загрузка выбранного файла
python scripts/encode_refs.py search --query 'fibroblast H3K27me3'
python scripts/encode_refs.py download --accession ENCFF... --label fibroblast_H3K27me3
```

### Тесты

```bash
python -m unittest discover tests
```

## Знаки и обозначения

- `logFC > 0` — сигнал выше в HS, `logFC < 0` — выше в BF.
- **HS_higher / BF_higher** — регионы, значимые по edgeR TREAT (FDR < 0.05, |log2FC| > 1).
- **BF_enriched / HS_enriched / common** — регионы пулированного MACS2 bdgdiff.
- **consensus** — объединение пулированных пиков BF и HS после фильтрации по blacklist. На нём считаются counts для edgeR.

## Данные и ссылки

- GEO: GSE40740.
- Pemberton H. et al. Genome-wide co-localization of Polycomb orthologs and their effects on gene expression in human fibroblasts. *Genome Biology* 15, R23 (2014) [doi:10.1186/gb-2014-15-2-r23](https://doi.org/10.1186/gb-2014-15-2-r23).
- Инструменты: MACS2, csaw, edgeR, ChIPseeker, deepTools, HOMER, IGV.

## Лицензия

MIT, см. [LICENSE](LICENSE).
