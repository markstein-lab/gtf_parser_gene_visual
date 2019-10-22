import matplotlib.pyplot as plt


class GeneImage(object):
    def __init__(self, gene_length, exon_intervals, exon_color="red", bar_color='gray', bg_color="white"):
        self.geneLength = gene_length
        self.ra_exons = exon_intervals[1]
        self.rb_exons = exon_intervals[2]

        self.exonColor = exon_color
        self.barColor= bar_color
        self.bgColor = bg_color

        self.ylims = {'exon_max': 1, 'exon_min':0.5} # ori = 2 and 1
        self.figure, self.canvas = plt.subplots(figsize=(13,4))
        
        self.canvas.set_facecolor(self.bgColor)
        self._draw()

    def _set_limits(self):
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


    def _draw(self):
        self._set_limits()

        # ra_exon
        for i in range(len(self.ra_exons)):
            self._draw_ra_exons(self.ra_exons[i][1])

        # rb_exon
        for i in range(len(self.rb_exons)):
            self._draw_rb_exons(self.rb_exons[i][1])

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
        
        self.canvas.text(1619000, 1.35, 'ABCB7', fontsize=12, ha='center')
        self.canvas.text(1619000, 1.9, 'ABCB7-RA', fontsize=12, ha='center')
        self.canvas.text(1619000, 2.9, 'ABCB7-RB', fontsize=12, ha='center')

        plt.axis('off')
        plt.savefig("output/output.png", bbox_inches='tight')
    
    def show(self):
        plt.show()
        