init:
	@python3 version.py
	@python3 -m venv env
	@./env/bin/pip3 install -r requirements.txt
	@mkdir -p input_files
	@mkdir -p bankit_files
	@mkdir -p study_scripts


freeze:
	@./env/bin/pip3 freeze > requirements.txt

ipython:
	@./env/bin/ipython

validate:
	@./env/bin/python3 validate_meta_data.py ./input_files/source_modifier.csv

align_seq:
	@./env/bin/python3 align_seq.py ./input_files/unaligned_sequences.fasta ./input_files/aligned_seq.json ./input_files/aligned_meta.csv

dump_seq_meta:
	@./env/bin/python3 dump_seq_meta.py ./input_files/aligned_seq.json ./input_files/aligned_meta.csv ./input_files/aligned_ignore.csv

generate: validate align_seq dump_seq_meta
	@./env/bin/python3 generate_bankit.py ./input_files/aligned_meta.csv ./input_files/source_modifier.csv ./bankit_files/bankit.fasta
	@./env/bin/python3 generate_feature.py ./input_files/aligned_meta.csv ./bankit_files/features.txt
	@./env/bin/python3 generate_sc.py ./input_files/source_modifier.csv ./input_files/treatment.csv ./bankit_files/structured_comments.csv

.PHONY:
	init freeze ipython validate align_seq dump_seq_meta generate

.DEFAULT_GOAL := generate
