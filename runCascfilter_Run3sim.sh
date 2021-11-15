#! /usr/bin/env bash

o2-analysis-timestamp -b --configuration json://${PWD}/triggerjson_Run3sim.json \
| o2-analysis-event-selection -b --configuration json://${PWD}/triggerjson_Run3sim.json \
| o2-analysis-multiplicity-table -b --configuration json://${PWD}/triggerjson_Run3sim.json \
| o2-analysis-trackextension -b --configuration json://${PWD}/triggerjson_Run3sim.json \
| o2-analysis-trackselection -b --configuration json://${PWD}/triggerjson_Run3sim.json \
| o2-analysis-centrality-table -b --configuration json://${PWD}/triggerjson_Run3sim.json \
| o2-analysis-pid-tpc -b --configuration json://${PWD}/triggerjson_Run3sim.json \
| o2-analysis-pid-tof -b --configuration json://${PWD}/triggerjson_Run3sim.json \
| o2-analysis-weak-decay-indices -b --configuration json://${PWD}/triggerjson_Run3sim.json \
| o2-analysis-lf-lambdakzerobuilder -b --configuration json://${PWD}/triggerjson_Run3sim.json \
| o2-analysis-lf-lambdakzeroanalysis -b --configuration json://${PWD}/triggerjson_Run3sim.json \
| o2-analysis-lf-cascadebuilder -b --configuration json://${PWD}/triggerjson_Run3sim.json \
| o2-analysis-lf-strangeness-filter -b --configuration json://${PWD}/triggerjson_Run3sim.json \
| o2-analysis-lf-cascadeanalysis -b --configuration json://${PWD}/triggerjson_Run3sim.json 