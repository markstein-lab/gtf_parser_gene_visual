import matplotlib.pyplot as plt


class GeneImage(object):
    def __init__(self, gene_name, gene_length, exon_intervals, ra_start_stop_codon, rb_start_stop_codon, exon_color="orange", bar_color='gray', bg_color="white", start_codon_color = "green", stop_codon_color = "red"):
        self.genename = gene_name.upper()

        self.geneLength = gene_length
        self.ra_exons = exon_intervals[1]
        self.rb_exons = exon_intervals[2]

        self.ra_start_stop_codon = ra_start_stop_codon
        self.rb_start_stop_codon = rb_start_stop_codon

        self.exonColor = exon_color
        self.barColor= bar_color
        self.bgColor = bg_color
        self.stopcodoncolor = stop_codon_color
        self.startcodoncolor = start_codon_color

        self.ylims = {'exon_max': 1, 'exon_min':0.5}
        self.figure, self.canvas = plt.subplots(figsize=(13,4))
        
        self.canvas.set_facecolor(self.bgColor)
        self._draw()


    def _set_bar_limits(self):
        self.ylims['bar_min'] = self.ylims['exon_max']+0.2
        self.ylims['bar_max'] = self.ylims['bar_min']+(self.ylims['exon_max']-self.ylims['exon_min'])/5.0
        

    # fill the canvas with ra
    def _draw_ra_exons(self, span):
        self.canvas.fill_between(span, self.ylims['exon_min']+1.5, self.ylims['exon_max']+1.5,
                                 edgecolor='black', facecolor=self.exonColor)
        return True


    # fill the canvas with rb
    def _draw_rb_exons(self, span):
        self.canvas.fill_between(span, self.ylims['exon_min']+2.5, self.ylims['exon_max']+2.5,
                                 edgecolor='black', facecolor=self.exonColor)
        return True


    # fill the canvas with RA start_codon and arrow
    def _draw_ra_start_codon(self, span):
        self.canvas.fill_between([span[0],span[0]], self.ylims['exon_min']+1.5, self.ylims['exon_max']+1.6, edgecolor='black', facecolor='black')
        self.canvas.arrow(span[0],2.6,100,0, head_width=0.1, head_length=20, linewidth=0.5, color='black')
        self.canvas.fill_between([span[1],span[1]], self.ylims['exon_min']+1.5, self.ylims['exon_max']+1.5, edgecolor=self.startcodoncolor, facecolor=self.startcodoncolor)
        return True
    

    # fill the canvas with RA stop_codon
    def _draw_ra_stop_codon(self, span):
        self.canvas.fill_between(span, self.ylims['exon_min']+1.5, self.ylims['exon_max']+1.5,
                                 edgecolor=self.stopcodoncolor, facecolor=self.stopcodoncolor)
        return True


    # fill the canvas with RB start_codon and arrow
    def _draw_rb_start_codon(self, span):
        self.canvas.fill_between([span[0],span[0]], self.ylims['exon_min']+2.5, self.ylims['exon_max']+2.6, edgecolor='black', facecolor='black')
        self.canvas.arrow(span[0],3.6,100,0, head_width=0.1, head_length=20, linewidth=0.5, color='black')
        self.canvas.fill_between([span[1],span[1]], self.ylims['exon_min']+2.5, self.ylims['exon_max']+2.5, edgecolor=self.startcodoncolor, facecolor=self.startcodoncolor)
        return True
    

    # fill the canvas with RB stop_codon
    def _draw_rb_stop_codon(self, span):
        self.canvas.fill_between(span, self.ylims['exon_min']+2.5, self.ylims['exon_max']+2.5,
                                 edgecolor=self.stopcodoncolor, facecolor=self.stopcodoncolor)
        return True


    def _draw(self):
        self._set_bar_limits()

        # ra_exon
        for i in range(len(self.ra_exons)):
            self._draw_ra_exons(self.ra_exons[i][1])

        # rb_exon
        for i in range(len(self.rb_exons)):
            self._draw_rb_exons(self.rb_exons[i][1])
        
        '''
        (Pdb) ra_start_stop_codon
        [('start_codon', [1618501, 1618503]), ('stop_codon', [1622987, 1622989])]
        (Pdb) rb_start_stop_codon 
        [('start_codon', [1618501, 1618503]), ('stop_codon', [1622987, 1622989])]
        '''
        # draw RA start_stop codon
        for i in self.ra_start_stop_codon:
            if i[0].endswith('start_codon'):
                self._draw_ra_start_codon(i[1])
            if i[0].endswith('stop_codon'):
                self._draw_ra_stop_codon(i[1])
        for j in self.rb_start_stop_codon:
            if j[0].endswith('start_codon'):
                self._draw_rb_start_codon(j[1])
            if j[0].endswith('stop_codon'):
                self._draw_rb_stop_codon(j[1])

        # draw gene line 
        self.canvas.fill_between([self.geneLength[0], self.geneLength[1]],
                                  self.ylims['bar_min'], self.ylims['bar_max'],
                                  edgecolor=self.bgColor, facecolor=self.barColor)
        # draw gene line for ra
        self.canvas.fill_between([self.geneLength[0], self.geneLength[1]],
                                  self.ylims['bar_min']+1, self.ylims['bar_max']+1,
                                  edgecolor=self.bgColor, facecolor=self.barColor)
        
        # draw gene line for rb
        self.canvas.fill_between([self.geneLength[0], self.geneLength[1]],
                                  self.ylims['bar_min']+2, self.ylims['bar_max']+2,
                                  edgecolor=self.bgColor, facecolor=self.barColor)

        self.canvas.text(self.geneLength[0], 1.35, self.genename, fontsize=12, ha='center')
        self.canvas.text(self.geneLength[0], 1.8, self.genename + '-RA', fontsize=12, ha='center')
        self.canvas.text(self.geneLength[0], 2.8, self.genename +'-RB', fontsize=12, ha='center')

        plt.axis('off')
        plt.savefig("output/test.png", bbox_inches='tight')
    

    def show(self):
        plt.show()
        