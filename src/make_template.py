import json
MODULE = 'deseq2'

mi_template_json = {'module_version': '00.00.02', 'program_name': 'deseq2', 'program_subname': '', 'program_version': '3.8', 'compute': {'environment': 'aws', 'language': 'Python', 'language_version': '3.7', 'vcpus': 2, 'memory': 6000}, 'program_arguments': '', 'program_input': [{'input_type': 'file', 'input_file_type': 'TXT', 'input_position': -1, 'input_prefix': '-i'}, {'input_type': 'file', 'input_file_type': 'CSV', 'input_position': -1, 'input_prefix': '-i'}], 'program_output': [{'output_type': 'folder', 'output_file_type': '', 'output_position': 0, 'output_prefix': '-o'}, {'output_type': 'file', 'output_file_type': 'CSV', 'output_position': -100, 'output_prefix': '-o'}], 'alternate_inputs': [{'input_type': 'file', 'input_file_type': 'CSV', 'input_position': 0, 'input_prefix': '-group'}], 'alternate_outputs': [], 'defaults': {"output_file": ""}}
with open(MODULE+'.template.json','w') as fout:
    json.dump(mi_template_json, fout)

io_dryrun_json = {'input': ['s3://hubseq-data/test/dnaseq_targeted/run_test1/expressionqc/expressionqc.counts_matrix.column.csv'], 'output': ['s3://hubseq-data/test/dnaseq_targeted/run_test1/deseq2/'],  'alternate_inputs': ['s3://hubseq-data/test/dnaseq_targeted/run_test1/expressionqc/expressionqc.samplegroups.csv'], 'alternate_outputs': [], 'program_arguments': '', 'sample_id': MODULE+'_test', 'dryrun': ''}
io_json = {'input': ['s3://hubseq-data/test/dnaseq_targeted/run_test1/expressionqc/expressionqc.counts_matrix.column.csv'], 'output': ['s3://hubseq-data/test/dnaseq_targeted/run_test1/deseq2/'],  'alternate_inputs': ['s3://hubseq-data/test/dnaseq_targeted/run_test1/expressionqc/expressionqc.samplegroups.csv'], 'alternate_outputs': [], 'program_arguments': '', 'sample_id': MODULE+'_test'}

with open(MODULE+'.dryrun_test.io.json','w') as fout:
    json.dump(io_dryrun_json, fout)
with open(MODULE+'.test.io.json','w') as fout:
    json.dump(io_json, fout)

# job info test JSONs                                                                                                        
job_json = {"container_overrides": {"command": ["--module_name", MODULE, "--run_arguments", "s3://hubseq-data/modules/"+MODULE+"/job/"+MODULE+".test.job.json", "--working_dir", "/home/"]}, "jobqueue": "batch_scratch_queue", "jobname": "job_"+MODULE+"_test"}
with open(MODULE+'.test.job.json','w') as fout:
    json.dump(io_json, fout)
