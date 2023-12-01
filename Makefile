init:
	@python -m venv env
	@./env/bin/pip3 install -r requirements.txt

freeze:
	@./env/bin/pip3 freeze > requirements.txt

ipython:
	@./env/bin/ipython

generate:
	@./env/bin/python3 generate.py

.PHONY:
	init freeze ipython

.DEFAULT_GOAL := generate
