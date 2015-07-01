#!/bin/bash
./bin/MINA -io sample/clinical.txt -ip sample/CNA.txt -ip sample/METH.txt -ip sample/mRNA.txt -s sample/symbol.txt -o output/ -m omi dist -dhi 0.3
