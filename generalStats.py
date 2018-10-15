import numpy as np
import plotly.graph_objs as go
import os, pickle

stats_dir = os.path.dirname(os.path.realpath(__file__)) + '/statistics'

mm10_genes = open(stats_dir + '/mm10_coding_genes.txt').readlines()
mm10_genes = [s.strip() for s in mm10_genes]
mm10_choices = [{'label': s, 'value': s} for s in mm10_genes]

hg38_genes = open(stats_dir + '/hg38_coding_genes.txt').readlines()
hg38_genes = [s.strip() for s in hg38_genes]
hg38_choices = [{'label': s, 'value': s} for s in hg38_genes]

##
# 
##
class GenomeStatistic:
  X = None
  Y = None
  unnormY = None
  step = None
  maxY = None
  xlabel = None

  def __init__(self, celltype, stat_nm, xlabel = 'stat_nm'):
    # Import statistics from file
    with open('%s/%s_%s.pkl' % (stats_dir, celltype, stat_nm), 'rb') as f:
      ans = pickle.load(f, encoding = 'latin1')
    BarBins, Ns, step = ans

    self.X = BarBins
    self.unnormY = Ns
    self.Y = Ns / sum(Ns)
    self.step = step
    self.maxY = max(self.Y)
    if xlabel == 'stat_nm':
      self.xlabel = stat_nm
    else:
      self.xlabel = xlabel
    return

  def cumulative(self, xval):
    for idx in range(len(self.X)):
      if self.X[idx] > xval:
        break
    if idx == 0:
      cum = 0
    else:
      cum = 100 * sum(self.Y[:idx-1])

    if cum <= 5:
      var_text, var_color = 'very low', 'rgb(221, 46, 31)'
    elif cum >= 95:
      var_text, var_color = 'very high', 'rgb(0, 140, 221)'
    elif cum <= 25:
      var_text, var_color = 'low', 'rgb(236, 100, 12)'
    elif cum >= 75:
      var_text, var_color = 'high', 'rgb(96, 170, 20)'
    else:
      var_text, var_color = 'typical', 'rgb(115, 118, 121)'
    return '%.1f' % (cum), var_text, var_color

  ## Barplot
  def trace(self, xval, color_light = 'rgb(158, 202, 225)', color_dark = 'rgb(108, 152, 175)'):
    colors = []
    for x in self.X:
      if x < xval:
        colors.append(color_dark)
      else:
        colors.append(color_light)
    return go.Bar(
      x = self.X,
      y = self.Y,
      text = ['Cumulative: %.2f%%' % (100*sum(self.Y[:idx])) for idx in range(len(self.Y))],
      opacity = 0.6,
      width = 0.01, # touching
      marker = dict(
        color = colors,
        line = dict(
          width = 0,
        ),
      ),
    )

  ##
  # Layout
  ##
  def layout(self, xval):
    ticks = 'outside'
    ticklen = 3
    tickwidth = 0.5
    max_x = max(xval, max(self.X))
    min_x = min(xval, min(self.X))
    return go.Layout(
      xaxis = dict(
        range = [min_x, max_x],
        title = self.xlabel,
        titlefont = dict(
          family = 'Arial',
        ),
        ticks = ticks,
        nticks = 10,
        ticklen = ticklen,
        tickwidth = tickwidth,
        fixedrange = True,
      ),
      yaxis = dict(
        autorange = True,
        showgrid = False,
        showline = False,
        ticks = '',
        showticklabels = False,
        fixedrange = True,
      ),
      font = dict(
        family = 'Arial',
      ),
      shapes = [
        dict(
          type = 'line',
          x0 = xval,
          y0 = 0,
          x1 = xval,
          y1 = self.maxY,
          line = dict(
            color = 'rgb(33, 33, 33)',
            width = 1.75,
          ),
        ),
      ],
      margin = dict(
        t = 20,
        # Right-aligned margin
        r = 10,
        l = 150,
        b = 70,
      ),
    )


# Tooltips
def get_tooltip_precision(var_text):
  if 'low' in var_text:
    return 'Low precision is highly stochastic. DNA repair outcomes are heterogeneous with no single repair genotype dominating.'
  elif 'high' in var_text:
    return 'High precision is like gambling with strongly loaded dice. DNA repair outcomes are homogeneous and enriched for just a handful of unique genotypes.'
  elif 'typical' in var_text:
    return 'Precision reflects the homogeneity or heterogeneity of predicted repair outcomes. High precision repair predominantly yields one or a handful of genotypes.'

def get_tooltip_phi(var_text):
  if 'low' in var_text:
    return 'Target sites with low microhomology strength reflect that most local microhomologies are short, GC-poor, and far from the cutsite. Repair outcomes will tend to be microhomology-less deletions and insertions.'
  elif 'high' in var_text:
    return 'Target sites with high microhomology strength have some local microhomologies that are long, GC-rich, or close to the cutsite. Repair outcomes will tend to be microhomology-based deletions.'
  elif 'typical' in var_text:
    return 'This microhomology strength score for the target site integrates all local microhomologies into a single number. Microhomology strength depends on distance to the cutsite, microhomology length and GC content.'

def get_tooltip_frameshift(var_text):
  if 'low' in var_text:
    return 'A low frameshift frequency will tend to keep a protein-coding gene in frame. The typical genomic frameshift frequency is above 66% because 1-bp insertions and 1-2 bp deletions are particularly common repair outcomes.'
  elif 'high' in var_text:
    return 'A high frameshift frequency will tend to knock a protein-coding gene out of frame. The typical genomic frameshift frequency is above 66% because 1-bp insertions and 1-2 bp deletions are particularly common repair outcomes.'
  elif 'typical' in var_text:
    return 'Frameshift frequency is informative for designing gRNAs for gene knockouts. The typical genomic frameshift frequency is above 66% because 1-bp insertions and 1-2 bp deletions are particularly common repair outcomes.'

##
# Init
##
celltypes = [
  'mESC',
  'U2OS',
  'HEK293',
  'HCT116',
  'K562',
]

stats = {
  'Precision': 'Precision score',
  'Phi': 'Microhomology strength score',
  'Frameshift frequency': 'Frameshift frequency',
}

GSD = dict()
for celltype in celltypes:
  for stat in stats:
    xlabel = stats[stat]
    gs = GenomeStatistic(celltype, stat, xlabel = xlabel)
    GSD[(celltype, stat)] = gs

# gs_mESC_precision = GenomeStatistic('Precision', xlabel = 'Precision score')
# gs_mESC_logphi = GenomeStatistic('Phi', xlabel = 'Log phi (Microhomology strength)')
# gs_mESC_frameshift = GenomeStatistic('Frameshift frequency')

# gs_onebpfq = GenomeStatistic('1-bp ins frequency', xlabel = '1-bp insertion frequency')
# gs_mhfq = GenomeStatistic('MH del frequency', xlabel = 'Microhomology deletion frequency')
# gs_mhlessfq= GenomeStatistic('MHless del frequency', xlabel = 'Microhomology-less deletion frequency')

# gs_frame0 = GenomeStatistic('Frame +0 frequency')
# gs_frame1 = GenomeStatistic('Frame +1 frequency')
# gs_frame2 = GenomeStatistic('Frame +2 frequency')

# gs_highest = GenomeStatistic('Highest outcome frequency')
# gs_highestdel = GenomeStatistic('Highest del frequency', xlabel = 'Highest deletion frequency')
# gs_highestins = GenomeStatistic('Highest ins frequency', xlabel = 'Highest insertion frequency')

# gs_expectedindellength = GenomeStatistic('Expected indel length')
