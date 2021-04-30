OUT_DIR = out

$(shell mkdir -p $(OUT_DIR))

#install workspace
install:
	sudo apt install python3-pip python3.8-venv
	python3 -m venv md
	. md/bin/activate; pip install -Ur requirements.txt
	echo \n
	echo "Sucessfully installed virtual environment. Please run 'source md/bin/activate' to activate the environment"

task3:
	python3 src/task3.py 10

# remove virtual environment
uninstall:
	rm -rf md

clean:
	rm -rf $(OUT_DIR)

