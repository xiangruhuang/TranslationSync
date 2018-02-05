graph_types=[1, 2, 3]
folder=./results

draw:
	matlab -nodesktop -r "draw_figure($(graph_types), '$(folder)'); exit;"
	scp $(folder)/Graph_*.eps xrhuang@narsil-13.cs.utexas.edu:./public_html/figures/TranslationSync/$(folder)/
