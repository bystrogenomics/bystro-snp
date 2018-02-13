# bystro-snp [![Build Status](https://travis-ci.org/akotlar/bystro-snp.svg?branch=master)](https://travis-ci.org/akotlar/bystro-snp)

## TL;DR

A really fast, simple [SNP](http://www.pnas.org/content/114/10/E1923) pre-processor and annotator. Millions of variants per minute.

```shell
go get github.com/akotlar/bystro-snp && go install $_;

pigz -d -c in.snp.gz | bystro-snp --minGq .95 | pigz -c - > output 2> log.txt
```

<br>

## Description



Performs several important functions:
1) Splits multiallelics
2) Performs QC on variants: checks whether allele is ACTG, +ACTG, or -Int
3) Filters samples based on genotype quality
4) Calculates whether site is transition, transversion, or neither
5) Processes all available samples
    - calculates homozygosity, heterozygosity, missingness
    - labels samples as homozygous, heterozygous, or missing

<br>

## Publication

bystro-snp is used to pre-proces SNP files for [Bystro](https://bystro.io) ([github](https://github.com/akotlar/bystro))

If you use bystro-snp please cite https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1387-3 

<br>

## Performance
Millions of variants/rows per minute. Performance is dependent on the # of samples.

<br>

## Installation

```shell
go get github.com/akotlar/bystro-snp && go install $_;
```

<br>

## Use

Via pipe:
```shell
pigz -d -c in.snp.gz | bystro-snp --minGq .95 | pigz -c - > out.gz
```

Via ```inPath``` argument:
```shell
bystro-snp --inPath in.snp --minGq .95 " > out
```

<br>

## Output
```tsv
chrom <String>   pos <Int>   type <String[SNP|DEL|INS|MULTIALLELIC]>    ref <String>    alt <String>    trTv <Int[0|1|2]>     heterozygotes <String>     heterozygosity <Float64>    homozygotes <String>     homozygosity <Float64>     missingGenos <String>    missingness <Float64>    sampleMaf <Float64>
```

<br>

## Optional arguments

```shell
--minGq <Float>
```

Minimum genotype quality to keep (0 - 1)

<br>

```shell
--inPath /path/to/uncompressedFile.vcf
```

An input file path, to an uncompressed VCF file. Defaults to `stdin`

<br>

```shell
--errPath /path/to/log.txt
```

Where to store log messages. Defaults to `STDERR`

<br>

```shell
--emptyField "!"
```

Which value to assign to missing data. Defaults to `!`

<br>

```shell
--fieldDelimiter ";"
```

Which delimiter to use when joining multiple values. Defaults to `;`
