#!/bin/bash

../bin/eggs -u 1 < sim.vcf.gz | zcat > test.out
../bin/eggs -p 0.5 < sim.vcf.gz | zcat >> test.out
../bin/eggs -s < sim.vcf.gz | zcat >> test.out
../bin/eggs -d 0.7,0.05 < sim.vcf.gz | zcat >> test.out
../bin/eggs -g 0.001 < sim.vcf.gz | zcat >> test.out
../bin/eggs -t < empirical.vcf.gz 2>> test.out
../bin/eggs -a < sim.vcf.gz | zcat >> test.out
../bin/eggs -x < sim.vcf.gz | zcat >> test.out
../bin/eggs -r empirical.vcf.gz < sim.vcf.gz | zcat >> test.out
../bin/eggs -b empirical.vcf.gz < sim.vcf.gz | zcat >> test.out
../bin/eggs -b 0.1,0.2 < sim.vcf.gz | zcat >> test.out
../bin/eggs -m empirical.vcf.gz < sim.vcf.gz | zcat >> test.out
../bin/eggs -b 0.1,0.2 -u 1 -p 0.5 -d 0.7,0.05 -s < sim.vcf.gz | zcat >> test.out

