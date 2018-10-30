import pickle, copy, os, datetime, subprocess, json
from collections import defaultdict
import random
import numpy as np
import pandas as pd
from scipy.stats import entropy
import time
from io import StringIO

import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_table_experiments as dt
from dash.dependencies import Input, Output, State
import plotly.graph_objs as go
import flask
import plotly
from flask_caching import Cache

import inDelphi
import generalStats
import lib, header

from indelphi_app import app

# from rq import Queue
# from worker import conn
# rq = Queue(connection = conn)

# init
inDelphi.init_model()
try:
  os.mkdir('user-csvs/')
except FileExistsError:
  pass
else:
  subprocess.check_output('rm -rf user-csvs/*', shell = True)

# Set up flask caching
CACHE_CONFIG = {
  'CACHE_TYPE': 'redis',
  'CACHE_REDIS_URL': os.environ.get('REDIS_URL', 'localhost:6379')
}
cache = Cache()
cache.init_app(app.server, config = CACHE_CONFIG)
cache_timeout = 300

# Remove these plotly modebar buttons to limit interactivity
modebarbuttons_2d = ['zoom2d', 'pan2d', 'select2d', 'lasso2d', 'zoomIn2d', 'zoomOut2d', 'autoScale2d', 'resetScale2d', 'hoverClosestCartesian', 'hoverCompareCartesian', 'toggleSpikelines']

# Random default, which is cached on filesystem
default_left_text = ''.join([random.choice(list('ACGT')) for s in range(60)])
default_right_text = ''.join([random.choice(list('ACGT')) for s in range(4)]) + 'GG' + ''.join([random.choice(list('ACGT')) for s in range(54)])
if os.path.isfile('single_default.pkl'):
  subprocess.check_output('rm -rf single_default.pkl', shell = True)

## Parameters

###################################################################
###################################################################
##
# App layout
##
layout = html.Div([

  ##
  # Hidden divs for light data storage
  ##
  html.Div(
    [
      html.Div(
        id = 'S_hidden-pred-signal',
        children = 'init'
      ),
      html.Div(
        id = 'S_hidden-chosen-celltype',
        children = 'mESC'
      ),
      html.Div(
        id = 'S_hidden-pred-df-summary-signal',
        children = 'init'
      ),
      html.Div(
        id = 'S_hidden-cache-dsb-left',
        children = '%s' % (time.time())
      ),
      html.Div(
        id = 'S_hidden-cache-dsb-right',
        children = '%s' % (time.time())
      ),
      html.Div(
        id = 'S_hidden-cache-pam-left',
        children = '%s' % (time.time())
      ),
      html.Div(
        id = 'S_hidden-cache-pam-right',
        children = '%s' % (time.time())
      ),
      html.Div(
        id = 'S_hidden-cache-revcomp',
        children = '%s' % (time.time())
      ),

      dcc.Location(
        id = 'S_url',
        refresh = False,
      ),
    ],
    style = dict(
      display = 'none',
    ),
  ),

  ##
  # Header
  ##
  html.Div(
    [
      ###################################################
      # Upper header
      ###################################################
      header.get_navigation_header('single'),

      ###################################################
      # Sequence boxes
      ###################################################
      html.Div([
        # Left box
        html.Div(
          [
            dcc.Input(
              id = 'S_textbox1', 
              size = 28,
              value = default_left_text,
              type = 'text',
              autofocus = True,
              style = dict(
                direction = 'rtl',
                fontFamily = 'monospace',
                fontSize = 16,
                float = 'right',
              ),
            )
          ],
          className = 'dna_textbox',
        ),

        # Right box
        html.Div(
          [
            dcc.Input(
              id = 'S_textbox2', 
              size = 28,
              # value = 'GTGTGGCTGAAGGCATAGTAATTCTGA',
              value = default_right_text,
              type = 'text',
              style = dict(
                textAlign = 'left',
                fontFamily = 'monospace',
                fontSize = 16,
                float = 'left',
              ),
            ),
            html.A('🔃 strand',
              id = 'S_button-revcomp',
              style = dict(
                fontSize = 16,
                textDecoration = 'none',
                verticalAlign = 'middle',
                float = 'right',
                position = 'relative',
                transform = 'translateY(27%)',
              ),
            ),
          ],
          className = 'dna_textbox',
        ),
      ], 
        style = dict(
          verticalAlign = 'center',
          whiteSpace = 'nowrap',
          overflowX = 'auto',
        ),
      ),

      ###################################################
      # Row of PAM arrows 
      ###################################################
      html.Div(
        [
          html.Div(
            [
              # (Left) placeholder
              html.Div([],
                style = dict(
                  display = 'table-cell',
                  width = '25%',
                ),
              ),

              ###################################################
              # (Center) DSB arrows: Simple
              ###################################################
              html.Div([
                  html.A('◄',
                    id = 'S_button-dsb-left',
                    style = dict(
                      textDecoration = 'none',
                      fontFamily = 'monospace',
                      color = 'rgb(0, 174, 179)',
                      fontSize = 20,
                    ),
                  ),
                  '\tCut site\t',
                  html.A('►',
                    id = 'S_button-dsb-right',
                    style = dict(
                      textDecoration = 'none',
                      fontFamily = 'monospace',
                      color = 'rgb(0, 174, 179)',
                      fontSize = 20,
                    ),
                  ),
                ],
                style = dict(
                  textAlign = 'center',
                  display = 'table-cell',
                  width = '50%',
                )
              ),

              ###################################################
              # (Right) PAM arrows
              ###################################################
              html.Div([
                html.A('◄',
                  id = 'S_button-pam-left',
                  style = dict(
                    textDecoration = 'none',
                    fontFamily = 'monospace',
                    fontSize = 20,
                    color = 'rgb(237, 71, 149)',
                  ),
                ),
                dcc.Input(
                  id = 'S_textbox_pam', 
                  size = 5,
                  value = 'NGG',
                  type = 'text',
                  autofocus = False,
                  style = dict(
                    fontFamily = 'monospace',
                    fontSize = 14,
                    textAlign = 'center',
                    height = '18px',
                    width = '70px',
                    marginLeft = '5px',
                    marginRight = '5px',
                    transform = 'translateY(-2px)',
                  ),
                ),
                html.A('►',
                  id = 'S_button-pam-right',
                  style = dict(
                    textDecoration = 'none',
                    fontFamily = 'monospace',
                    fontSize = 20,
                    color = 'rgb(237, 71, 149)',
                  ),
                ),
              ],
              style = dict(
                display = 'table-cell',
                textAlign = 'right',
                width = '25%',
              )),
            ],
            style = dict(
              display = 'table-row',
            ),
          ),
        ], 
        style = dict(
          display = 'table',
          width = '100%',
          # Defaults to 32px tall for some reason, force it to 24 px
          marginTop = '-4px',
          marginBottom = '-4px',
        ),
      ),

      ################################################### 
      # Row for cell-type choice
      ###################################################
      html.Div(
        [
          html.Div(
            [
              # Left blank
              html.Div([],
                style = dict(
                  textAlign = 'center',
                  display = 'table-cell',
                  width = '35%',
                )
              ),

              # Middle
              html.Div(
                [
                  html.Span(
                    'Cell type: mESC',
                    id = 'S_celltype_chosen',
                  ),
                ],
                style = dict(
                  textAlign = 'center',
                  display = 'table-cell',
                  width = '30%',
                )
              ),

              # Right
              html.Div(
                [
                  html.A(
                    'HCT116',
                    id = 'S_celltype_link1',
                    n_clicks_timestamp = '0',
                    style = dict(color = 'gray'),
                  ),
                  ', ',
                  html.A(
                    'HEK293',
                    id = 'S_celltype_link2',
                    n_clicks_timestamp = '0',
                    style = dict(color = 'gray'),
                  ),
                  ', ',
                  html.A(
                    'K562',
                    id = 'S_celltype_link3',
                    n_clicks_timestamp = '0',
                    style = dict(color = 'gray'),
                  ),
                  ', ',
                  html.A(
                    'U2OS',
                    id = 'S_celltype_link4',
                    n_clicks_timestamp = '0',
                    style = dict(color = 'gray'),
                  ),
                ],
                style = dict(
                  textAlign = 'right',
                  display = 'table-cell',
                  width = '35%',
                )
              )
            ],
            style = dict(
              display = 'table-row',
            ),
          )
        ],
        style = dict(
          display = 'table',
          width = '100%',
        )
      )

    ],
    style = dict(
      position = 'fixed',
      top = 0,
      backgroundColor = 'white',
      borderBottom = '3px solid #777777',
      zIndex = 1e6,
      width = '1010px',
      left = '50%',
      transform = 'translate(-50%, 0px)',
    ),
  ),

  ##
  # Body / plots
  ##
  html.Div(
    [
      ###################################################
      # Module: Summary, alignment preview
      ###################################################
      html.Div([
        # header
        html.Div([
          html.Div([
            html.Strong(
              '',
              id = 'S_text-summary-module-header'
            )],
            className = 'module_header_text'),
          ],
          className = 'module_header'
        ),

        html.Div(
          [
            # Text table
            dcc.Graph(
              id = 'S_summary-alignment-table',
              config = dict(
                modeBarButtonsToRemove = modebarbuttons_2d,
                displaylogo = False,
                displayModeBar = False,
              ),
              style = dict(
                height = 290,
                width = 629,
              ),
              className = 'eight columns',
            ),

            # bar chart
            dcc.Graph(
              id = 'S_summary-alignment-barchart',
              config = dict(
                modeBarButtonsToRemove = modebarbuttons_2d,
                displaylogo = False,
                displayModeBar = False,
              ),
              style = dict(
                height = 300,
                width = 300,
              ),
              className = 'four columns',
            ),
          ],
          className = 'row',
        ),
      ], className = 'module_style',
      ),

      # Animate bottom
      html.Div([

        ###################################################
        # Module: Indel Lengths
        ###################################################
        html.Div([
          # header
          html.Div([
            html.Div([
              html.Strong('Indel length predictions')
              ],
              className = 'module_header_text'),
            ],
            className = 'module_header'
          ),

          html.Div(
            [
              # Frameshift
              html.Div(
                [
                  dcc.Graph(
                    id = 'S_plot-fs',
                    style = dict(
                      height = 300, 
                      width = 215,
                    ),
                    config = dict(
                      modeBarButtonsToRemove = modebarbuttons_2d,
                      displaylogo = False,
                      displayModeBar = False,
                    ),
                  ),
                ],
                className = 'three columns',
              ),

              # Indel length
              html.Div(
                [
                  dcc.Graph(
                    id = 'S_plot-indel-len',
                    style = dict(
                      height = 300, 
                      width = 715,
                    ),
                    config = dict(
                      modeBarButtonsToRemove = modebarbuttons_2d,
                      displaylogo = False,
                      displayModeBar = False,
                    ),
                  ),
                ],
                className = 'nine columns',
              ),
            ],
            className = 'row',
          ),

          html.Div([
            html.Div(
              html.A(
                '🔗 Shareable link to your results', 
                id = 'S_page-link'
              ),
            ),
          ], style = dict(
              height = '30px',
              textAlign = 'center',
            ),
          ),

        ], className = 'module_style',
        ),

        ###################################################
        # Module: Genome statistics
        ###################################################
        ## Precision
        html.Div([
          # header
          html.Div([
            html.Div([
              html.Strong('Comparison to predictions at {:,} SpCas9 target sites in human exons and introns'.format(int(sum(generalStats.GSD[('mESC', 'Phi')].unnormY))))
              ],
              className = 'module_header_text'),
            ],
            className = 'module_header'
          ),

          ## Precision
          html.Div(
            [
              html.Div(
                [
                  dcc.Graph(
                    id = 'S_plot-genstats-precision',
                    style = dict(
                      height = 200, 
                      width = 970/2,
                    ),
                    config = dict(
                      modeBarButtonsToRemove = modebarbuttons_2d,
                      displaylogo = False,
                      displayModeBar = False,
                    ),
                  ),
                ],
                className = 'six columns',
              ),
              html.Div(
                [
                  html.Div(
                    id = 'S_text-genstats-precision',
                    className = 'generalstats_text_inner',
                  ),
                ],
                style = dict(
                  height = 200,
                ),
                className = 'six columns',
              ),
            ],
            className = 'row',
          ),

          ## Phi
          html.Div(
            [
              html.Div(
                [
                  dcc.Graph(
                    id = 'S_plot-genstats-logphi',
                    style = dict(
                      height = 200, 
                      width = 970/2,
                    ),
                    config = dict(
                      modeBarButtonsToRemove = modebarbuttons_2d,
                      displaylogo = False,
                      displayModeBar = False,
                    ),
                  ),
                ],
                className = 'six columns',
              ),
              html.Div(
                [
                  html.Div(
                    id = 'S_text-genstats-logphi',
                    className = 'generalstats_text_inner',
                  ),
                ],
                style = dict(
                  height = 200,
                ),
                className = 'six columns',
              ),
            ],
            className = 'row',
          ),

          ## Frameshift frequency
          html.Div(
            [
              html.Div(
                [
                  dcc.Graph(
                    id = 'S_plot-genstats-frameshift',
                    style = dict(
                      height = 200, 
                      width = 970/2,
                    ),
                    config = dict(
                      modeBarButtonsToRemove = modebarbuttons_2d,
                      displaylogo = False,
                      displayModeBar = False,
                    ),
                  ),
                ],
                className = 'six columns',
              ),
              html.Div(
                [
                  html.Div(
                    id = 'S_text-genstats-frameshift',
                    className = 'generalstats_text_inner',
                  ),
                ],
                style = dict(
                  height = 200,
                ),
                className = 'six columns',
              ),
            ],
            className = 'row',
          ),

          html.Div([
            html.Div(
              html.A(
                '📑 Download summary statistics',
                id = 'S_summary-download-link'
              ),
            ),
          ], style = dict(
              height = '30px',
              textAlign = 'center',
            ),
          ),

        ], className = 'module_style',
        ),

        ###################################################
        # Module: New version of last module
        ###################################################
        html.Div([
          # header
          html.Div([
            html.Div([
              html.Strong('All predictions of 1-bp insertion and 1- to 60-bp deletion events')
              ],
              className = 'module_header_text'),
            ],
            className = 'module_header'
          ),

          # Settings
          # Row: Display indel type...
          html.Div(
            [
              html.Strong(
                'Display indel type:',
                style = dict(
                  textAlign = 'right',
                  marginRight = '5px',
                  height = '36px',  # height of one dropdown line
                  lineHeight = '36px',  # centers vertically
                ),
                className = 'three columns',
              ),

              # Multi drop down to select columns
              dcc.Dropdown(
                id = 'S_dropdown-indel-type',
                options = [
                  {'label': '1-bp insertions', 'value': '1-bp insertions'},
                  {'label': 'Microhomology deletions', 'value': 'Microhomology deletions'},
                  {'label': 'Microhomology-less deletions', 'value': 'Microhomology-less deletions'},
                ],
                multi = True,
                searchable = False,
                clearable = False,
                value = ['1-bp insertions', 'Microhomology deletions', 'Microhomology-less deletions'],
                className = 'nine columns',
              ),
            ],
            style = dict(
              marginBottom = '5px',
              marginTop = '10px',
            ),
            className = 'row',
          ),

          # Row: Sort by...
          html.Div(
            [
              html.Strong(
                'Sort by:',
                style = dict(
                  textAlign = 'right',
                  marginRight = '5px',
                  height = '36px',  # height of one dropdown line
                  lineHeight = '36px',  # centers vertically
                ),
                className = 'three columns',
              ),
              # Sorting columns
              dcc.Dropdown(
                id = 'S_dropdown-sort',
                options = [
                  {'label': 'Indel length', 'value': 'Indel length'},
                  {'label': 'Predicted frequency', 'value': 'Predicted frequency'},
                ],
                value = 'Indel length',
                searchable = False,
                clearable = False,
                className = 'four columns',
              ),
              html.Div('',
                className = 'five columns',
              ),
            ],
            style = dict(
              marginBottom = '5px',
              marginTop = '10px',
            ),
            className = 'row',
          ),

          # Row: Display frequencies within...
          html.Div(
            [
              html.Strong(
                'Frequency range:',
                style = dict(
                  textAlign = 'right',
                  marginRight = '5px',
                  height = '36px',  # height of one dropdown line
                  lineHeight = '36px',  # centers vertically
                ),
                className = 'three columns',
              ),
              # Sorting columns
              dcc.RangeSlider(
                id = 'S_rangeslider-freq',
                min = 0,
                max = 7,
                step = 1,
                value = [1, 7],
                marks = {
                  0: '0.05%',
                  1: '0.5%',
                  2: '1%',
                  3: '5%',
                  4: '10%',
                  5: '25%',
                  6: '50%',
                  7: '100%',
                },
                allowCross = False,
                className = 'eight columns',
              ),
              html.Div('',
                className = 'one columns',
              ),
            ],
            style = dict(
              marginBottom = '5px',
              marginTop = '10px',
            ),
            className = 'row',
          ),

          # Row: Display indel range...
          html.Div(
            [
              html.Strong(
                'Indel length range:',
                style = dict(
                  textAlign = 'right',
                  marginRight = '5px',
                  height = '36px',  # height of one dropdown line
                  lineHeight = '36px',  # centers vertically
                ),
                className = 'three columns',
              ),
              # Sorting columns
              dcc.RangeSlider(
                id = 'S_rangeslider-indel-len',
                min = 0,
                max = 36,
                step = 1,
                dots = True,
                value = [0, 36],
                marks = {i: '{} bp'.format(1 - i) if bool(1 - i <= 0) else '+1 bp' for i in [1 - 1, 1 - -1, 1 - -5, 1 - -10, 1 - -15, 1 - -20, 1 - -25, 1 - -30, 1 - -35]},
                allowCross = False,
                className = 'eight columns',
              ),
              html.Div('',
                className = 'one columns',
              ),
            ],
            style = dict(
              marginBottom = '5px',
              marginTop = '10px',
            ),
            className = 'row',
          ),

          # Download link
          html.Div([
            html.Div(
              html.A(
                '📜 Download table of inDelphi genotype predictions', 
                id = 'S_csv-download-link'
              ),
            ),
          ], style = dict(
              height = '30px',
              textAlign = 'center',
            ),
          ),

          # Table and barplot
          dcc.Graph(
            id = 'S_plot-table-genotypes-v2',
            config = dict(
              modeBarButtonsToRemove = modebarbuttons_2d,
              displaylogo = False,
              displayModeBar = False,
            ),
          ),

          html.Div(
            dcc.Markdown(
              'Copyright MIT 2018.\nAll Rights Reserved. [Terms of use](https://www.crisprindelphi.design/termsofuse)',
            ),
            style = dict(
              textAlign = 'center',
              height = '40px',
            )
          ),

        ], className = 'module_style',
        ),

        ],
        id = 'S_plots_body',
        style = dict(
          display = 'none',
        ),
        className = 'animate-bottom',
      ),

    ],
    # body style
    # id = 'S_plots_body',
    style = dict(
      # display = 'none',
      transform = 'translateY(%spx)' % (200),
    ),
  ),
  ##

  ],  # body div
  style = dict(
    width = '970px',
    margin = '0 auto',
  )
)

#######################################################################
#########################      CALLBACKS      #########################
#######################################################################

## Arrow buttons
# These must be 1->1 mapping, otherwise we can't tell which n_clicks triggered the callback
@app.callback(
  Output('S_hidden-cache-dsb-left', 'children'),
  [Input('S_button-dsb-left', 'n_clicks')])
def update_cache_dsb_left(n_clicks):
  return '%s' % (time.time())

@app.callback(
  Output('S_hidden-cache-dsb-right', 'children'),
  [Input('S_button-dsb-right', 'n_clicks')])
def update_cache_dsb_right(n_clicks):
  return '%s' % (time.time())

@app.callback(
  Output('S_hidden-cache-pam-left', 'children'),
  [Input('S_button-pam-left', 'n_clicks')])
def update_cache_pam_left(n_clicks):
  return '%s' % (time.time())

@app.callback(
  Output('S_hidden-cache-pam-right', 'children'),
  [Input('S_button-pam-right', 'n_clicks')])
def update_cache_pam_right(n_clicks):
  return '%s' % (time.time())

@app.callback(
  Output('S_hidden-cache-revcomp', 'children'),
  [Input('S_button-revcomp', 'n_clicks')])
def update_cache_revcomp(n_clicks):
  return '%s' % (time.time())

@app.callback(
  Output('S_textbox1', 'value'),
  [Input('S_hidden-cache-dsb-left', 'children'),
   Input('S_hidden-cache-dsb-right', 'children'),
   Input('S_hidden-cache-pam-left', 'children'),
   Input('S_hidden-cache-pam-right', 'children'),
   Input('S_hidden-cache-revcomp', 'children'),
   Input('S_url', 'pathname')],
  [State('S_textbox1', 'value'),
   State('S_textbox2', 'value'),
   State('S_textbox_pam', 'value')])
def update_textbox1_arrow(cache_dsb_left, cache_dsb_right, cache_pam_left, cache_pam_right, cache_rc, url, text1, text2, text_pam):
  left_dsb_time = float(cache_dsb_left)
  right_dsb_time = float(cache_dsb_right)
  left_pam_time = float(cache_pam_left)
  right_pam_time = float(cache_pam_right)
  rc_time = float(cache_rc)

  pageload_time_threshold = 0.03
  pageload_dsb = bool(abs(left_dsb_time - right_dsb_time) < pageload_time_threshold)
  pageload_pam = bool(abs(left_pam_time - right_pam_time) < pageload_time_threshold)
  pageload_rc = bool(abs(rc_time - left_pam_time) < pageload_time_threshold)
  if pageload_dsb and pageload_pam and pageload_rc:
    valid_flag, celltype, seq, cutsite = lib.parse_valid_url_path_single(url)
    if not valid_flag or cutsite is None:
      return text1
    else:
      return seq[:cutsite]

  latest_dsb_click = max(left_dsb_time, right_dsb_time)
  latest_pam_click = max(left_pam_time, right_pam_time)
  latest_rc_click = rc_time

  if latest_dsb_click > max(latest_pam_click, latest_rc_click):
    # Clicked DSB
    if left_dsb_time > right_dsb_time:
      return text1[:-1]
    elif right_dsb_time > left_dsb_time:
      return text1 + text2[0]
  elif latest_pam_click > max(latest_dsb_click, latest_rc_click):
    # Clicked PAM
    if left_pam_time > right_pam_time:
      newtext1, newtext2 = lib.pam_shift(text1, text2, text_pam, 'left')
      return newtext1
    elif right_pam_time > left_pam_time:
      newtext1, newtext2 = lib.pam_shift(text1, text2, text_pam, 'right')
      return newtext1
  elif latest_rc_click > max(latest_dsb_click, latest_pam_click):
    # Clicked RC button
    return lib.revcomp(text2)

@app.callback(
  Output('S_textbox2', 'value'),
  [Input('S_hidden-cache-dsb-left', 'children'),
   Input('S_hidden-cache-dsb-right', 'children'),
   Input('S_hidden-cache-pam-left', 'children'),
   Input('S_hidden-cache-pam-right', 'children'),
   Input('S_hidden-cache-revcomp', 'children'),
   Input('S_url', 'pathname')],
  [State('S_textbox1', 'value'),
   State('S_textbox2', 'value'),
   State('S_textbox_pam', 'value')])
def update_textbox2_arrow(cache_dsb_left, cache_dsb_right, cache_pam_left, cache_pam_right, cache_rc, url, text1, text2, text_pam):
  left_dsb_time = float(cache_dsb_left)
  right_dsb_time = float(cache_dsb_right)
  left_pam_time = float(cache_pam_left)
  right_pam_time = float(cache_pam_right)
  rc_time = float(cache_rc)
  
  pageload_time_threshold = 0.03
  pageload_dsb = bool(abs(left_dsb_time - right_dsb_time) < pageload_time_threshold)
  pageload_pam = bool(abs(left_pam_time - right_pam_time) < pageload_time_threshold)
  pageload_rc = bool(abs(rc_time - left_pam_time) < pageload_time_threshold)
  if pageload_dsb and pageload_pam and pageload_rc:
    valid_flag, celltype, seq, cutsite = lib.parse_valid_url_path_single(url)
    if not valid_flag or cutsite is None:
      return text2
    else:
      return seq[cutsite:]

  latest_dsb_click = max(left_dsb_time, right_dsb_time)
  latest_pam_click = max(left_pam_time, right_pam_time)
  latest_rc_click = rc_time

  if latest_dsb_click > max(latest_pam_click, latest_rc_click):
    # Clicked DSB
    if left_dsb_time > right_dsb_time:
      return text1[-1] + text2
    elif right_dsb_time > left_dsb_time:
      return text2[1:]
  elif latest_pam_click > max(latest_dsb_click, latest_rc_click):
    # Clicked PAM
    if left_pam_time > right_pam_time:
      newtext1, newtext2 = lib.pam_shift(text1, text2, text_pam, 'left')
      return newtext2
    elif right_pam_time > left_pam_time:
      newtext1, newtext2 = lib.pam_shift(text1, text2, text_pam, 'right')
      return newtext2
  elif latest_rc_click > max(latest_dsb_click, latest_pam_click):
    # Clicked RC button
    return lib.revcomp(text1)

##
# Celltype choice callbacks
## 
@app.callback(
  Output('S_hidden-chosen-celltype', 'children'),
  [Input('S_celltype_link1', 'n_clicks_timestamp'),
   Input('S_celltype_link2', 'n_clicks_timestamp'),
   Input('S_celltype_link3', 'n_clicks_timestamp'),
   Input('S_celltype_link4', 'n_clicks_timestamp'),
   Input('S_url', 'pathname')],
  [State('S_hidden-chosen-celltype', 'children')])
def update_hidden_celltype(link1, link2, link3, link4, url, prev_state):
  pageload_time_threshold = 0.03
  links = [int(s) for s in [link1, link2, link3, link4]]
  if abs(max(links) - min(links)) <= pageload_time_threshold:
    valid_flag, celltype, seq, cutsite = lib.parse_valid_url_path_single(url)
    if not valid_flag:
      return prev_state
    else:
      return celltype

  celltypes = ['HCT116', 'HEK293', 'K562', 'mESC', 'U2OS']
  celltypes.remove(prev_state)
  clicked_idx = links.index(max(links))
  return celltypes[clicked_idx]

@app.callback(
  Output('S_celltype_chosen', 'children'),
  [Input('S_hidden-chosen-celltype', 'children')])
def update_celltype_chosen_text(celltype):
  return 'Cell type: %s' % (celltype)

@app.callback(
  Output('S_celltype_link1', 'children'),
  [Input('S_hidden-chosen-celltype', 'children')])
def update_celltype_link1(celltype):
  celltypes = ['HCT116', 'HEK293', 'K562', 'mESC', 'U2OS']
  celltypes.remove(celltype)
  return '%s' % (celltypes[0])

@app.callback(
  Output('S_celltype_link2', 'children'),
  [Input('S_hidden-chosen-celltype', 'children')])
def update_celltype_link2(celltype):
  celltypes = ['HCT116', 'HEK293', 'K562', 'mESC', 'U2OS']
  celltypes.remove(celltype)
  return '%s' % (celltypes[1])

@app.callback(
  Output('S_celltype_link3', 'children'),
  [Input('S_hidden-chosen-celltype', 'children')])
def update_celltype_link3(celltype):
  celltypes = ['HCT116', 'HEK293', 'K562', 'mESC', 'U2OS']
  celltypes.remove(celltype)
  return '%s' % (celltypes[2])

@app.callback(
  Output('S_celltype_link4', 'children'),
  [Input('S_hidden-chosen-celltype', 'children')])
def update_celltype_link4(celltype):
  celltypes = ['HCT116', 'HEK293', 'K562', 'mESC', 'U2OS']
  celltypes.remove(celltype)
  return '%s' % (celltypes[3])

##
# Module header callbacks
##
@app.callback(
  Output('S_text-summary-module-header', 'children'),
  [Input('S_textbox1', 'value'),
   Input('S_textbox2', 'value')])
def update_summary_module_header(text1, text2):
  presumed_grna = text1[-17:] + text2[:3]
  return 'Summary of predictions at target site with gRNA: %s' % (presumed_grna)

##
# Prediction callback
##
@cache.memoize(timeout = cache_timeout)
def indelphi_predict_cache(seq, cutsite, celltype):
  # if default
  if seq == default_left_text + default_right_text and cutsite == 60 and celltype == 'mESC':
    pkl_fn = 'single_default.pkl'
    if os.path.isfile(pkl_fn):
      ans = pickle.load(open(pkl_fn, 'rb'))
      pred_df, stats = ans
    else:
      pred_df, stats = inDelphi.predict(seq, int(cutsite), celltype)
      stats = pd.DataFrame(stats, index = [0])
      with open(pkl_fn, 'wb') as f:
        pickle.dump(tuple([pred_df, stats]), f)
  # If not default
  else:
    pred_df, stats = inDelphi.predict(seq, int(cutsite), celltype)
    stats = pd.DataFrame(stats, index = [0])
  return pred_df, stats

@app.callback(
  Output('S_hidden-pred-signal', 'children'),
  [Input('S_textbox1', 'value'),
   Input('S_textbox2', 'value'),
   Input('S_hidden-chosen-celltype', 'children')])
def update_pred_df(text1, text2, celltype):
  seq = text1 + text2
  seq = seq.upper()
  cutsite = len(text1)
  indelphi_predict_cache(seq, cutsite, celltype)
  return '%s,%s,%s' % (seq, cutsite, celltype)

## 
# Style callbacks: Hide figures on page load until processing is complete
##
@app.callback(
  Output('S_plots_body', 'style'),
  [Input('S_plot-genstats-precision', 'figure')],
  # [Input('S_summary-alignment-barchart', 'figure')],
  [State('S_plots_body', 'style')])
def update_plots_body_style(fig, prev_style):
  new_style = prev_style
  if fig is not None:
    if 'display' in new_style:
      del new_style['display']
  return new_style

##
# Summary of predictions callbacks
##
def indelphi_summary_cache(signal):
  seq, cutsite, celltype = signal.split(',')
  pred_df, stats = indelphi_predict_cache(seq, cutsite, celltype)

  mhless_gt_df = inDelphi.add_mhless_genotypes(pred_df, stats, length_cutoff = 6)
  inDelphi.add_genotype_column(mhless_gt_df, stats)
  top10 = mhless_gt_df.sort_values('Predicted frequency', ascending = False).iloc[:10]
  return top10

@app.callback(
  Output('S_hidden-pred-df-summary-signal', 'children'),
  [Input('S_hidden-pred-signal', 'children')])
def update_pred_df_top10_summary(signal):
  indelphi_summary_cache(signal)
  return signal

# change
@app.callback(
  Output('S_summary-alignment-table', 'figure'),
  [Input('S_hidden-pred-df-summary-signal', 'children'),
   Input('S_hidden-pred-signal', 'children'),
  ])
def update_summary_alignment_text(summary_signal, signal):
  top10 = indelphi_summary_cache(summary_signal)

  seq, cutsite, celltype = signal.split(',')
  pred_df, stats = indelphi_predict_cache(seq, cutsite, celltype)
  
  # mhless_gt_df = inDelphi.add_mhless_genotypes(pred_df, stats)
  # top10 = mhless_gt_df.sort_values('Predicted frequency', ascending = False).iloc[:10]

  # top10 = pred_df.sort_values('Predicted frequency', ascending = False).iloc[:10]
  # inDelphi.add_genotype_column(pred_df, stats)
  # inDelphi.add_genotype_column(mhless_gt_df, stats)

  gts = top10['Genotype']
  fqs = top10['Predicted frequency']
  lens = top10['Length']
  cats = top10['Category']
  fq_strings = ['-'] + ['%.1f' % (s) for s in fqs]

  gap_gts = lib.get_gapped_alignments(top10, stats)
  
  alignments, categories = [], []
  cutsite = stats['Cutsite'].iloc[0]
  reference_seq = stats['Reference sequence'].iloc[0]
  reference_seq = lib.add_bar(reference_seq, cutsite)
  alignments.append(lib.trim_alignment(reference_seq, cutsite, 'ref'))
  categories.append('Reference')

  for gt, length, cat in zip(gap_gts, lens, cats):
    alignments.append(gt)
    if cat == 'ins':
      categories.append('%s-bp insertion' % (length))
    elif cat == 'del':
      categories.append('%s-bp deletion' % (length))

  return dict(
    data = [go.Table(
      columnwidth = [310, 110, 60],
      header = dict(
        values = ['Alignment', 'Category', '%'],
        align = ['center', 'right', 'right'], 
        line = dict(color = 'white'),
        fill = dict(color = 'white'),
      ),
      cells = dict(
        values = [alignments, categories, fq_strings],
        align = ['center', 'right', 'right'], 
        line = dict(color = 'white'),
        fill = dict(color = 'white'),
        font = dict(family = 'monospace'),
      )
    )],
    layout = go.Layout(
      font = dict(
        family = 'monospace',
      ),
      margin = dict(
        l = 10,
        r = 0,
        t = 5,
        b = 5,
      ),
    ),
  )

# change
@app.callback(
  Output('S_summary-alignment-barchart', 'figure'),
  [Input('S_hidden-pred-df-summary-signal', 'children'),
   Input('S_hidden-pred-signal', 'children'),
  ])
def update_summary_alignment_barchart(summary_signal, signal):
  top10 = indelphi_summary_cache(summary_signal)
  seq, cutsite, celltype = signal.split(',')
  pred_df, stats = indelphi_predict_cache(seq, cutsite, celltype)

  # mhless_gt_df = inDelphi.add_mhless_genotypes(pred_df, stats)
  # inDelphi.add_genotype_column(mhless_gt_df, stats)
  # top10 = mhless_gt_df.sort_values('Predicted frequency', ascending = False).iloc[:10]

  # inDelphi.add_genotype_column(pred_df, stats)
  # top10 = pred_df.sort_values('Predicted frequency', ascending = False).iloc[:10]
  fqs = top10['Predicted frequency'][::-1]

  cats = top10['Category'][::-1]
  mhls = top10['Microhomology length'][::-1]
  colors = []
  for cat, mhl in zip(cats, mhls):
    if mhl == 0:
      colors.append('rgb(236, 100, 12)')
    else:
      if cat == 'del':
        colors.append('rgb(221, 46, 31)')
      else:
        colors.append('rgb(0, 190, 220)')

  return dict(
    data = [go.Bar(
      x = fqs,
      y = np.arange(len(fqs)),
      orientation = 'h',
      marker = dict(
        color = colors,
        line = dict(
          width = 0,
        ),
      ),
    )],
    layout = go.Layout(
      xaxis = dict(
        fixedrange = True,
      ),
      yaxis = dict(
        showticklabels = False,
        fixedrange = True,
      ),
      margin = dict(
        l = 5,
        r = 10,
        t = 58,
        b = 42,
      ),
    ),
  )

##
# General stats callbacks
##
@app.callback(
  Output('S_plot-genstats-precision', 'figure'),
  [Input('S_hidden-pred-signal', 'children')])
def plot_genstats_precision(signal):
  seq, cutsite, celltype = signal.split(',')
  pred_df, stats = indelphi_predict_cache(seq, cutsite, celltype)

  xval = stats['Precision'].iloc[0]
  return dict(
    data = [
      generalStats.GSD[(celltype, 'Precision')].trace(xval),
    ],
    layout = generalStats.GSD[(celltype, 'Precision')].layout(xval),
  )

@app.callback(
  Output('S_plot-genstats-logphi', 'figure'),
  [Input('S_hidden-pred-signal', 'children')])
def plot_genstats_logphi(signal):
  seq, cutsite, celltype = signal.split(',')
  pred_df, stats = indelphi_predict_cache(seq, cutsite, celltype)

  xval = np.log(stats['Phi'].iloc[0])
  return dict(
    data = [
      generalStats.GSD[(celltype, 'Phi')].trace(xval),
    ],
    layout = generalStats.GSD[(celltype, 'Phi')].layout(xval),
  )

@app.callback(
  Output('S_plot-genstats-frameshift', 'figure'),
  [Input('S_hidden-pred-signal', 'children')])
def plot_genstats_frameshift(signal):
  seq, cutsite, celltype = signal.split(',')
  pred_df, stats = indelphi_predict_cache(seq, cutsite, celltype)

  xval = stats['Frameshift frequency'].iloc[0]
  return dict(
    data = [
      generalStats.GSD[(celltype, 'Frameshift frequency')].trace(xval),
    ],
    layout = generalStats.GSD[(celltype, 'Frameshift frequency')].layout(xval),
  )


## General stats text
@app.callback(
  Output('S_text-genstats-precision', 'children'),
  [Input('S_hidden-pred-signal', 'children')])
def text_genstats_precision(signal):
  seq, cutsite, celltype = signal.split(',')
  pred_df, stats = indelphi_predict_cache(seq, cutsite, celltype)

  xval = stats['Precision'].iloc[0]
  cum, var_text, var_color = generalStats.GSD[(celltype, 'Precision')].cumulative(xval)
  tooltip_msg = generalStats.get_tooltip_precision(var_text)
  return [
    html.Strong('This target site has '),
    html.Strong(var_text,
      style = dict(color = var_color),
    ),
    html.Strong(' precision. ',
      style = dict(color = var_color),
    ),
    html.Div(
      [
        html.Img(
          src = '/staticfiles/tooltip_logo',
          className = 'tooltiplogo',
        ),
        html.Span(
          tooltip_msg,
          className = 'tooltiptext'
        ),
      ], 
      className = 'tooltip',
    ),
    html.Br(),
    html.Span('Precision score: %.2f' % (xval),
      className = 'generalstats_subtext_style'),
    html.Br(),
    html.Span('Percentile: %s' % (cum),
      className = 'generalstats_subtext_style'),
  ]

@app.callback(
  Output('S_text-genstats-logphi', 'children'),
  [Input('S_hidden-pred-signal', 'children')])
def text_genstats_logphi(signal):
  seq, cutsite, celltype = signal.split(',')
  pred_df, stats = indelphi_predict_cache(seq, cutsite, celltype)

  xval = np.log(stats['Phi'].iloc[0])
  cum, var_text, var_color = generalStats.GSD[(celltype, 'Phi')].cumulative(xval)
  tooltip_msg = generalStats.get_tooltip_phi(var_text)
  return [
    html.Strong('This target site has '),
    html.Strong(var_text,
      style = dict(color = var_color),
    ),
    html.Strong(' microhomology strength. ',
      style = dict(color = var_color),
    ),
    html.Div(
      [
        html.Img(
          src = '/staticfiles/tooltip_logo',
          className = 'tooltiplogo',
        ),
        html.Span(
          tooltip_msg,
          className = 'tooltiptext'
        ),
      ], 
      className = 'tooltip',
    ),
    html.Br(),
    html.Span('Microhomology strength score: %.2f' % (xval),
      className = 'generalstats_subtext_style'),
    html.Br(),
    html.Span('Percentile: %s' % (cum),
      className = 'generalstats_subtext_style'),
  ]

@app.callback(
  Output('S_text-genstats-frameshift', 'children'),
  [Input('S_hidden-pred-signal', 'children')])
def text_genstats_frameshift(signal):
  seq, cutsite, celltype = signal.split(',')
  pred_df, stats = indelphi_predict_cache(seq, cutsite, celltype)

  xval = stats['Frameshift frequency'].iloc[0]
  cum, var_text, var_color = generalStats.GSD[(celltype, 'Frameshift frequency')].cumulative(xval)
  tooltip_msg = generalStats.get_tooltip_frameshift(var_text)
  return [
    html.Strong('This target site has '),
    html.Strong(var_text,
      style = dict(color = var_color),
    ),
    html.Strong(' frameshift frequency. ',
      style = dict(color = var_color),
    ),
    html.Div(
      [
        html.Img(
          src = '/staticfiles/tooltip_logo',
          className = 'tooltiplogo',
        ),
        html.Span(
          tooltip_msg,
          className = 'tooltiptext'
        ),
      ], 
      className = 'tooltip',
    ),
    html.Br(),
    html.Span('Frameshift frequency: %.1f%%' % (xval),
      className = 'generalstats_subtext_style'),
    html.Br(),
    html.Span('Percentile: %s' % (cum),
      className = 'generalstats_subtext_style'),
  ]

##
# Indel length and frameshift callbacks
@app.callback(
  Output('S_plot-indel-len', 'figure'),
  [Input('S_hidden-pred-signal', 'children'),
  ])
def plot_indel_len(signal):
  seq, cutsite, celltype = signal.split(',')
  pred_df, stats = indelphi_predict_cache(seq, cutsite, celltype)

  if False:
    lendf = inDelphi.get_indel_length_fqs(pred_df)
    X, Y = [], []
    for idx, row in lendf.iterrows():
      idl = int(row['Indel length'])
      if idl >= -35:
        X.append(idl)
        Y.append(row['Predicted frequency'])
    
    barmode_flag = None
    fig_data = [
      go.Bar(
        x = X,
        y = Y,
        opacity = 0.6,
        width = 0.9, 
        marker = dict(
          color = 'rgb(200, 20, 20)',
          line = dict(
            width = 0,
          ),
        ),
      )
    ]
  elif True:
    bdf = inDelphi.get_indel_length_breakdown(pred_df)
    X, Y = defaultdict(list), defaultdict(list)
    for indel_len in ['+1'] + [str(s) for s in range(-1, -35 - 1, -1)]:
      idl = str(indel_len)
      for detail in ['A', 'C', 'G', 'T', 'MH-less', 'Microhomology']:
        crit = (bdf['Detail'] == detail) & (bdf['Indel length'] == idl)
        ss = bdf[crit]
        if len(ss) == 0:
          yval = 0
        else:
          yval = sum(bdf[crit]['Predicted frequency'])
        if 'M' not in detail and indel_len != '+1':
          continue
        if 'M' in detail and indel_len == '+1':
          continue
        X[detail].append(indel_len)
        Y[detail].append(yval)

    colors = {
      'A': '#7CB82F',
      'C': '#00AEB3',
      'G': '#68C7EC',
      'T': '#00A0DC',
      'Microhomology': 'rgb(221, 46, 31)',
      'MH-less': 'rgb(236, 100, 12)',
    }

    traces = []
    for detail in ['A', 'C', 'G', 'T', 'MH-less', 'Microhomology']:
      trace = go.Bar(
        x = X[detail],
        y = Y[detail],
        name = detail,
        opacity = 0.6,
        width = 0.9,
        marker = dict(
          color = colors[detail],
          line = dict(width = 0),
        ),
      )
      traces.append(trace)
    
    barmode_flag = 'stack'
    fig_data = traces

  return dict(
    data = fig_data,
    layout = go.Layout(
      barmode = barmode_flag,
      showlegend = False,
      xaxis = dict(
        autorange = 'reversed',
        title = 'Indel length',
        ticks = 'outside',
        ticklen = 3,
        tickwidth = 0.5,
        tick0 = 1,
        dtick = 2,
        zeroline = False,
        fixedrange = True,
      ),
      yaxis = dict(
        title = 'Frequency (%)',
        hoverformat = '.2f%%',
        zeroline = False,
        fixedrange = True,
      ),
      font = dict(
        family = 'Arial',
      ),
      margin = dict(
        t = 30,
        l = 70,
        r = 20,
        b = 70,
      ),
    ),
  )

@app.callback(
  Output('S_plot-fs', 'figure'),
  [Input('S_hidden-pred-signal', 'children'),
  ])
def plot_fs(signal):
  seq, cutsite, celltype = signal.split(',')
  pred_df, stats = indelphi_predict_cache(seq, cutsite, celltype)

  fs_df = inDelphi.get_frameshift_fqs(pred_df)
  X = ['+0', '+1', '+2']
  Y = [float(fs_df[fs_df['Frame'] == s]['Predicted frequency']) for s in X]
  return dict(
    data = [
      go.Bar(
        x = X,
        y = Y,
        text = ['%.0f%%' % (s) for s in Y],
        textfont = dict(
          family = 'Arial',
        ),
        textposition = 'auto',
        opacity = 0.6,
        # width = 1,  # touching
        marker = dict(
          color = 'rgb(158, 202, 225)',
          line = dict(
            color = 'rgb(8, 48, 107)',
            width = 1.5,
          ),
        ),
      ),
    ],
    layout = go.Layout(
      # title = 'Frameshift frequency (%)',
      # margin = go.Margin(l = 40, r = 0, t = 40, b = 30),
      xaxis = dict(
        title = 'Frame',
        showline = False,
        tickvals = list(range(3)),
        ticktext = ['+0', '+1', '+2'],
        fixedrange = True,
      ),
      yaxis = dict(
        title = 'Frameshift frequency (%)',
        titlefont = dict(
          family = 'Arial',
        ),
        zeroline = False,
        showline = True,
        ticks = 'outside',
        tick0 = 0,
        ticklen = 3,
        tickwidth = 0.5,
        hoverformat = '.2f%%',
        fixedrange = True,
      ),
      font = dict(
        family = 'Arial',
      ),
      margin = dict(
        t = 30,
        l = 70,
        r = 20,
        b = 70,
      ),
    ),
  )

##
# Genotype table v2 callbacks
##
@app.callback(
  Output('S_plot-table-genotypes-v2', 'figure'),
  [Input('S_hidden-pred-signal', 'children'),
   Input('S_dropdown-indel-type', 'value'),
   Input('S_dropdown-sort', 'value'),
   Input('S_rangeslider-freq', 'value'),
   Input('S_rangeslider-indel-len', 'value'),
  ])
def update_genotype_table_v2(signal, indel_types, sort_col, freq_range, indel_len_range):
  seq, cutsite, celltype = signal.split(',')
  pred_df, stats = indelphi_predict_cache(seq, cutsite, celltype)

  # Filter indel types and by indel length range
  filter_1bp_ins = bool(indel_len_range[0] > 0)
  if '1-bp insertions' not in indel_types or filter_1bp_ins:
    pred_df = pred_df[pred_df['Category'] != 'ins']
  ins_crit = (pred_df['Category'] == 'ins')
  if 'Microhomology deletions' not in indel_types:
    crit = (pred_df['Category'] == 'del') & (pred_df['Microhomology length'] == 0)
    pred_df = pred_df[ins_crit | crit]
  if 'Microhomology-less deletions' not in indel_types:
    crit = (pred_df['Category'] == 'del') & (pred_df['Microhomology length'] > 0)
    pred_df = pred_df[ins_crit | crit]

  # Filter del len range
  [min_dellen, max_dellen] = [max(s - 1, 0) for s in indel_len_range]
  crit = (pred_df['Category'] == 'del') & (pred_df['Length'] >= min_dellen) & (pred_df['Length'] <= max_dellen)
  pred_df = pred_df[ins_crit | crit]

  # Expand MHless genotypes
  mhless_gt_df = inDelphi.add_mhless_genotypes(pred_df, stats, length_cutoff = indel_len_range[1])

  # Filter frequency range
  freq_mapper = {
    0: 0.05,
    1: 0.5,
    2: 1,
    3: 5,
    4: 10,
    5: 25,
    6: 50,
    7: 100,
  }
  [min_freq, max_freq] = [freq_mapper[s] for s in freq_range]
  crit = (mhless_gt_df['Predicted frequency'] >= min_freq) & (mhless_gt_df['Predicted frequency'] <= max_freq)
  mhless_gt_df = mhless_gt_df[crit]

  # Sort and process into table
  if sort_col == 'Predicted frequency':
    mhless_gt_df = mhless_gt_df.sort_values('Predicted frequency', ascending = True)
  else:
    tdf2 = mhless_gt_df[mhless_gt_df['Category'] == 'del'].sort_values('Length', ascending = False)
    temp_df = mhless_gt_df[mhless_gt_df['Category'] == 'ins']
    temp_df = tdf2.append(temp_df, ignore_index = True)
    mhless_gt_df = temp_df
  mhless_gt_df.reset_index(inplace = True)
  inDelphi.add_genotype_column(mhless_gt_df, stats)

  ins_colors = {
    'A': '#7CB82F',
    'C': '#00AEB3',
    'G': '#68C7EC',
    'T': '#00A0DC',
  }

  gap_gts = lib.get_gapped_alignments(mhless_gt_df, stats)

  cutsite = stats['Cutsite'].iloc[0]
  reference_seq = stats['Reference sequence'].iloc[0]
  reference_seq = lib.add_bar(reference_seq, cutsite)
  reference_seq = lib.trim_alignment(reference_seq, cutsite, 'ref')
  reference_ticktext = '%s  Indel len.  MH len.  Freq.   ' % (reference_seq)

  # Build bar plot with alignment text values for ticks 
  X, yticktexts, colors = [], [], []
  for idx, row in mhless_gt_df.iterrows():
    # Add color
    if row['Category'] == 'del':
      if row['Microhomology length'] == 0:
        colors.append('rgb(236, 100, 12)')
      else:
        colors.append('rgb(221, 46, 31)')
      indel_len = '-%s bp' % (int(row['Length']))
      mhbp = row['Microhomology length']
      if np.isnan(mhbp):
        mhbp = 0
      else:
        mhbp = '%s nt' % (int(mhbp))

    else:
      colors.append(ins_colors[row['Inserted Bases']])
      indel_len = '+1 bp'
      mhbp = '----'

    # Add x, y
    freq = row['Predicted frequency']
    X.append(freq)

    fixedwidth_mhbp = mhbp
    while len(fixedwidth_mhbp) != 5:
      fixedwidth_mhbp = ' ' + fixedwidth_mhbp

    fixedwidth_freq = '%.2f%%' % (freq)
    while len(fixedwidth_freq) != 6:
      fixedwidth_freq = ' ' + fixedwidth_freq

    fixedwidth_indel_len = indel_len
    while len(fixedwidth_indel_len) != 6:
      fixedwidth_indel_len = ' ' + fixedwidth_indel_len

    yticktexts.append('%s    %s    %s   %s   ' % (gap_gts[idx], fixedwidth_indel_len, fixedwidth_mhbp, fixedwidth_freq))

  ref_trace = go.Bar(
    x = [0],
    y = [1],
    orientation = 'h',
    opacity = 0.6,
    hoverinfo = 'x',
    width = 0.9,
    marker = dict(
      color = colors,
      line = dict(width = 0),
    ),
  )

  main_trace = go.Bar(
    x = X,
    y = np.arange(len(X)),
    orientation = 'h',
    opacity = 0.6,
    hoverinfo = 'x',
    width = 0.9,
    marker = dict(
      color = colors,
      line = dict(width = 0),
    ),
  )

  fig = plotly.tools.make_subplots(rows = 2, cols = 1)
  fig.append_trace(ref_trace, 1, 1)
  fig.append_trace(main_trace, 2, 1)

  fig['layout']['showlegend'] = False
  fig['layout']['hovermode'] = 'y'
  fig['layout']['hoverlabel'].update(
    font = dict(
      size = 12,
      family = 'monospace',
    )
  )
  fig['layout']['xaxis1'].update(
    showline = False,
    zeroline = False,
    showgrid = False,
    showticklabels = False,
    hoverformat = '.2f%%',
    fixedrange = True,
  )
  fig['layout']['yaxis1'].update(
    tickvals = [1],
    ticktext = [reference_ticktext],
    tickfont = dict(
      size = 12,
      family = 'monospace'
    ),
    domain = [1 - (1 / len(mhless_gt_df)), 1],
    fixedrange = True,
  )
  fig['layout']['xaxis2'].update(
    title = 'Frequency (%)',
    hoverformat = '.2f%%',
    zeroline = False,
    fixedrange = True,
  )
  fig['layout']['yaxis2'].update(
    tickvals = np.arange(len(yticktexts)),
    ticktext = yticktexts,
    tickfont = dict(
      size = 12,
      family = 'monospace'
    ),
    fixedrange = True,
    domain = [0, 1 - 1.5 * (1 / len(mhless_gt_df))],
  )
  fig['layout']['height'] = max(105 + 14 + 14 * len(mhless_gt_df), 200)
  fig['layout']['width'] = 920
  fig['layout']['margin'] = {
    'l': 600,
    'r': 25,
    't': 25,
    'b': 80,
  }
  return fig

##
# Download callbacks
##
@app.callback(
  Output('S_csv-download-link', 'href'), 
  [Input('S_hidden-pred-signal', 'children')])
def update_link(signal):
  seq, cutsite, celltype = signal.split(',')
  pred_df, stats = indelphi_predict_cache(seq, cutsite, celltype)

  pdf = inDelphi.add_mhless_genotypes(pred_df, stats)
  # inDelphi.add_genotype_column(pred_df, stats)
  # inDelphi.add_name_column(pred_df, stats)
  inDelphi.add_genotype_column(pdf, stats)
  inDelphi.add_name_column(pdf, stats)

  time = str(datetime.datetime.now()).replace(' ', '_').replace(':', '-')
  link_fn = '/dash/urlToDownload?value={}'.format(time)
  # pred_df.to_csv('user-csvs/%s.csv' % (time))
  pdf.to_csv('user-csvs/%s.csv' % (time))
  return link_fn

@app.callback(
  Output('S_summary-download-link', 'href'),
  [Input('S_hidden-pred-signal', 'children')],
  [State('S_page-link', 'href')])
def update_summary_link(signal, pagelink):
  seq, cutsite, celltype = signal.split(',')
  pred_df, stats = indelphi_predict_cache(seq, cutsite, celltype)

  stats['URL'] = pagelink

  time = str(datetime.datetime.now()).replace(' ', '_').replace(':', '-')
  link_fn = '/dash/urlToDownloadSummary?value={}'.format(time)
  stats.to_csv('user-csvs/%s_summary.csv' % (time))
  return link_fn

##
# Flask serving
##
@app.server.route('/dash/urlToDownload') 
def download_csv():
  value = flask.request.args.get('value')
  # create a dynamic csv or file here using `StringIO` 
  # (instead of writing to the file system)
  local_csv_fn = value.split('/')[-1]
  return flask.send_file(
    open('user-csvs/%s.csv' % (local_csv_fn), 'rb'),
    mimetype = 'text/csv',
    attachment_filename = 'inDelphi_output.csv',
    as_attachment = True,
  )

@app.server.route('/dash/urlToDownloadSummary') 
def download_summary_csv():
  value = flask.request.args.get('value')
  # create a dynamic csv or file here using `StringIO` 
  # (instead of writing to the file system)
  local_csv_fn = value.split('/')[-1]
  return flask.send_file(
    open('user-csvs/%s_summary.csv' % (local_csv_fn), 'rb'),
    mimetype = 'text/csv',
    attachment_filename = 'inDelphi_targetsite_summary.csv',
    as_attachment = True,
  )

@app.server.route('/staticfiles/tooltip_logo')
def serve_image():
  # BE VERY CAREFUL NOT TO SERVE ARBITRARY FILES
  return flask.send_from_directory(os.getcwd() + '/staticfiles/', 'noun_646495_cc.png')


##
# Page link callback
##
@app.callback(
  Output('S_page-link', 'href'),
  [Input('S_textbox1', 'value'),
   Input('S_textbox2', 'value'),
   Input('S_hidden-chosen-celltype', 'children')
  ])
def update_pagelink(text1, text2, celltype):
  seq = text1 + text2
  cutsite = len(text1)
  return 'https://www.crisprindelphi.design%s' % (lib.encode_dna_to_url_path_single(seq, cutsite, celltype))

