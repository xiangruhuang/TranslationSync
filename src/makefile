all:
	g++ -o solver main.cpp -O3 -std=c++11 -fopenmp
	g++ -o dataI artsquad.cpp -O3 -std=c++11 -fopenmp

n=2000
edge_density=0.2
bias=1.0
inc=0
max_iter=10000
noise_type=2
noise_ratio=0
decay=0.9
sigma=0.01
output=uniform_n1p2_$(sigma)/
a=-1
b=2
stopping=0.1
graph_file=

synthetic:
	mkdir -p $(output)
	./solver $(n) $(edge_density) $(bias) $(inc) $(noise_type) $(noise_ratio) $(max_iter) $(decay) 	$(a)	$(b)	$(sigma)	$(output)		$(graph_file)

graph1:
	mkdir -p $(output)/$@/
	./solver 2000 	0.1 			1.0 	0.0 	2 		 $(noise_ratio) 	10000 		0.8 	-1		2		$(sigma)	$(output)/$@/	graph1.meta ${stopping}

graph2:
	mkdir -p $(output)/$@/
	./solver 2000 	0.4 			0.2 	0.6 	2 		$(noise_ratio) 		10000 		0.8 	-1		2		$(sigma)	$(output)/$@/	graph2.meta ${stopping}

graph3:
	mkdir -p $(output)/$@/
	./solver 20000 	0.003 			1.0 	0.0 	2 			$(noise_ratio) 	10000 		0.8 	-1		2		$(sigma)	$(output)/$@/   graph3.meta ${stopping}

graph4:
	mkdir -p $(output)/$@/
	./solver 20000 	0.1 			0.07 	0.21 	2 				$(noise_ratio) 	10000 		0.8 	-1		2		$(sigma)	$(output)/$@/ 	graph4.meta ${stopping}

draw_graph:
	python graph.py $(sigma) 1 graph
	python graph.py $(sigma) 2 graph
	#python graph.py $(sigma) 3 graph
	python graph.py $(sigma) 4 graph

draw_time:
	python graph.py $(sigma) 1 time
	python graph.py $(sigma) 2 time
	#python graph.py $(sigma) 3 time
	python graph.py $(sigma) 4 time

