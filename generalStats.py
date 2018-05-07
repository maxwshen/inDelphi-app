import numpy as np
import plotly.graph_objs as go
import os, pickle

stats_dir = os.path.dirname(os.path.realpath(__file__)) + '/statistics'

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

  def __init__(self, stat_nm, xlabel = 'stat_nm'):
    # Import statistics from file
    with open('%s/%s.pkl' % (stats_dir, stat_nm), 'rb') as f:
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
        ticklen = ticklen,
        tickwidth = tickwidth,
      ),
      yaxis = dict(
        autorange = True,
        showgrid = False,
        showline = False,
        ticks = '',
        showticklabels = False,
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
        t = 10,
      ),
    )

##
# Init
##
gs_precision = GenomeStatistic('Precision', xlabel = 'Precision score')
gs_logphi = GenomeStatistic('logphi', xlabel = 'Log phi (Microhomology strength)')

gs_onebpfq = GenomeStatistic('1-bp ins frequency', xlabel = '1-bp insertion frequency')
gs_mhfq = GenomeStatistic('MH del frequency', xlabel = 'Microhomology deletion frequency')
gs_mhlessfq= GenomeStatistic('MHless del frequency', xlabel = 'Microhomology-less deletion frequency')

gs_frameshift = GenomeStatistic('Frameshift frequency')
gs_frame0 = GenomeStatistic('Frame +0 frequency')
gs_frame1 = GenomeStatistic('Frame +1 frequency')
gs_frame2 = GenomeStatistic('Frame +2 frequency')

gs_highest = GenomeStatistic('Highest outcome frequency')
gs_highestdel = GenomeStatistic('Highest del frequency', xlabel = 'Highest deletion frequency')
gs_highestins = GenomeStatistic('Highest ins frequency', xlabel = 'Highest insertion frequency')

gs_expectedindellength = GenomeStatistic('Expected indel length')
