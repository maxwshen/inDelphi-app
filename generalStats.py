import numpy as np
import plotly.graph_objs as go

trace_precision = None

##
# IO
##
def import_statistics(stat_nm):
  # Load from file
  X = np.arange(0, 1, .01)
  Y = [1/s for s in range(1, 101)]
  Y = np.array(Y) / sum(Y)
  return X, Y

##
# Generic template
def trace_template(stats, color):
  X, Y = stats
  return go.Bar(
    x = X,
    y = Y,
    text = ['Cumulative: %.2f%%' % (100*sum(Y[:idx])) for idx in range(len(Y))],
    opacity = 0.6,
    width = 0.01, # touching
    marker = dict(
      color = color,
      line = dict(
        width = 0,
      ),
    ),
  )

##
# Layout
##
def layout(xlabel, xval):
  ticks = 'outside'
  ticklen = 3
  tickwidth = 0.5
  return go.Layout(
    xaxis = dict(
      range = [0, 1],
      title = xlabel,
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
      showline = True,
      ticks = ticks,
      ticklen = ticklen,
      tickwidth = tickwidth,
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
        y1 = 0.1,
        line = dict(
          color = 'rgb(0, 0, 0)',
          width = 1.5,
        ),
      ),
    ],
  )

##
# Init
##
def init_all_traces():
  # Set global variables

  # Traces
  global trace_precision
  stats = import_statistics('precision')
  trace_precision = trace_template(stats, 'rgb(158, 202, 225)')

  # etc

  return