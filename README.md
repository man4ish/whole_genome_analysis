
# Whole Genome

## Validating wdl file
java -jar bin/womtool-59.jar validate whole_genome_analysis.wdl

## Generating input.json for workflow
java -jar bin/womtool-59.jar inputs whole_genome_analysis.wdl > input.json

## Running wdl locally on cromwell engine
java -jar bin/cromwell-58.jar run whole_genome_analysis.wdl --inputs input.json


# Related Links 
### Reference for WDL docs: 
### https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md
### https://dockstore.org/
### https://docs.docker.com/docker-hub/
### https://github.com/gatk-workflows/
### https://github.com/openwdl/learn-wdl

