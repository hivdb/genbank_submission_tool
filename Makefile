init:
	@python3 version.py
	@python -m venv env
	@./env/bin/pip3 install -r requirements.txt
	@mkdir -p input_files
	@mkdir -p bankit_files


freeze:
	@./env/bin/pip3 freeze > requirements.txt

ipython:
	@./env/bin/ipython

align_seq:
	@./env/bin/python3 align_seq.py ./input_files/unaligned_sequences.fasta ./input_files/aligned_seq.json ./input_files/aligned_meta.csv

generate: align_seq
	@./env/bin/python3 generate.py ./input_files/aligned_sequences.json ./input_files/source_modifier.csv

.PHONY:
	init freeze ipython

.DEFAULT_GOAL := generate
