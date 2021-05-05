OUT_DIR = out

$(shell mkdir -p $(OUT_DIR))


all: task2 task3 task4

#install workspace
install:
	sudo apt install python3-pip python3.8-venv
	python3 -m venv md
	. md/bin/activate; pip install -Ur requirements.txt
	echo \n
	echo "Sucessfully installed virtual environment. Please run 'source md/bin/activate' to activate the environment"

task2:
	python3 src/task2.py 200 6.3 1 out/task2.txt

task3:
	python3 src/task3.py out/task2.txt 0.001 10

task4:
	python3 src/task4.py out/task3.txt

# remove virtual environment
uninstall:
	rm -rf md

clean:
	rm -rf $(OUT_DIR)

