init:
	@python3 version.py
	@python -m venv env
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

generate: validate align_seq
	@./env/bin/python3 generate.py ./input_files/aligned_meta.csv ./input_files/source_modifier.csv "Human immunodeficiency virus 1" "Homo sapiens" ./bankit_files/bankit.fasta ./bankit_files/features.txt

.PHONY:
	init freeze ipython

.DEFAULT_GOAL := generate
