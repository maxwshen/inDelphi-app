import pickle, copy, os, datetime, subprocess
from collections import defaultdict
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

import inDelphi
import generalStats
import lib, header

from indelphi_app import app

# init
inDelphi.init_model()
if not os.path.isdir('user-csvs/'):
  os.mkdir('user-csvs/')
else:
  subprocess.check_output('rm -rf user-csvs/*', shell = True)

# Remove these plotly modebar buttons to limit interactivity
modebarbuttons_2d = ['zoom2d', 'pan2d', 'select2d', 'lasso2d', 'zoomIn2d', 'zoomOut2d', 'autoScale2d', 'resetScale2d', 'hoverClosestCartesian', 'hoverCompareCartesian', 'toggleSpikelines']

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
        id = 'S_hidden-pred-df',
        children = 'init'
      ),
      html.Div(
        id = 'S_hidden-pred-stats',
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
      # html.H4(
      #   'inDelphi',
      #   style = dict(
      #     textAlign = 'center',
      #   ),
      # ),
      header.navigation_header,

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
              value = 'TAGTTTCTAGCACAGCGCTGGTGTGGC',
              type = 'text',
              autofocus = True,
              style = dict(
                textAlign = 'right',
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
              value = 'GTGTGGCTGAAGGCATAGTAATTCTGA',
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
              html.Div('',
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
                  '\tDSB\t',
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
                    verticalAlign = 'middle',
                  ),
                ),
                dcc.Input(
                  id = 'S_textbox_pam', 
                  size = 5,
                  value = 'NGG',
                  type = 'text',
                  autofocus = True,
                  style = dict(
                    fontFamily = 'monospace',
                    fontSize = 14,
                    textAlign = 'center',
                    height = '18px',
                    width = '70px',
                    marginLeft = '5px',
                    marginRight = '5px',
                  ),
                ),
                html.A('►',
                  id = 'S_button-pam-right',
                  style = dict(
                    textDecoration = 'none',
                    fontFamily = 'monospace',
                    fontSize = 20,
                    color = 'rgb(237, 71, 149)',
                    verticalAlign = 'middle',
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
        ),
      ),

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
              html.Strong('Comparison to predictions at {:,} SpCas9 target sites in human exons and introns'.format(int(sum(generalStats.gs_logphi.unnormY))))
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
        ], className = 'module_style',
        ),

        ###################################################
        # Module: Detailed genotypes
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

          dcc.Graph(
            id = 'S_plot-table-genotypes',
            style = dict(
              width = 970,
            ),
            config = dict(
              modeBarButtonsToRemove = modebarbuttons_2d,
              displaylogo = False,
              displayModeBar = False,
            ),
          ),

          dt.DataTable(
            id = 'S_table-genotypes',
            rows = [{}], # init rows
            row_selectable = True,
            filterable = True,
            sortable = True,
            editable = False,
            selected_row_indices = [],
            column_widths = [970/4] * 4,
          ),

          html.Div([
            html.Div(
              html.A(
                '📑 Download target site summary and statistics',
                id = 'S_summary-download-link'
              ),
            ),
            html.Div(
              html.A(
                '📜 Download table of inDelphi genotype predictions', 
                id = 'S_csv-download-link'
              ),
            ),
            html.Div(
              html.A(
                '🔗 Shareable link to your results', 
                id = 'S_page-link'
              ),
            ),
          ], style = dict(
              marginLeft = '20',
            ),
          ),

          html.Div(
            'Copyright MIT 2018.\nAll Rights Reserved.',
            style = dict(
              textAlign = 'center',
              marginTop = '30',
              marginBottom = '30',
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
      transform = 'translateY(%spx)' % (160),
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
    valid_flag, seq, cutsite = lib.parse_valid_url_path_single(url)
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
    valid_flag, seq, cutsite = lib.parse_valid_url_path_single(url)
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
@app.callback(
  Output('S_hidden-pred-df', 'children'),
  [Input('S_textbox1', 'value'),
   Input('S_textbox2', 'value')])
def update_pred_df(text1, text2):
  seq = text1 + text2
  seq = seq.upper()
  cutsite = len(text1)
  pred_df, stats = inDelphi.predict(seq, cutsite)
  return pred_df.to_csv()

@app.callback(
  Output('S_hidden-pred-stats', 'children'),
  [Input('S_textbox1', 'value'),
   Input('S_textbox2', 'value')])
def update_pred_stats(text1, text2):
  seq = text1 + text2
  cutsite = len(text1)
  pred_df, stats = inDelphi.predict(seq, cutsite)
  return pd.DataFrame(stats, index = [0]).to_csv()

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
@app.callback(
  Output('S_summary-alignment-table', 'figure'),
  [Input('S_hidden-pred-df', 'children'),
   Input('S_hidden-pred-stats', 'children'),
  ])
def update_summary_alignment_text(pred_df_string, pred_stats_string):
  pred_df = pd.read_csv(StringIO(pred_df_string), index_col = 0)
  stats = pd.read_csv(StringIO(pred_stats_string), index_col = 0)
  
  inDelphi.add_genotype_column(pred_df, stats)

  top10 = pred_df.sort_values('Predicted frequency', ascending = False).iloc[:10]
  gts = top10['Genotype']
  fqs = top10['Predicted frequency']
  lens = top10['Length']
  cats = top10['Category']
  fq_strings = ['-'] + ['%.1f' % (s) for s in fqs]

  def trim_alignment(gt, cutsite, name):
    radius = 26
    if name == 'ins':
      trim_cand = gt[cutsite - radius : cutsite + radius + 1]
      if len(trim_cand) == 2*radius + 1:
        return trim_cand
      else:
        return gt
    else:
      trim_cand = gt[cutsite - radius : cutsite + radius + 1]
      if len(trim_cand) == 2*radius + 1:
        return trim_cand
      else:
        return gt
    return

  def add_bar(seq, cutsite):
    return seq[:cutsite] + '|' + seq[cutsite:]

  def get_gapped_alignments(top, stats):
    cutsite = stats['Cutsite'].iloc[0]
    gapped_aligns = []
    for idx, row in top.iterrows():
      gt = row['Genotype']
      gt_pos = row['Genotype position']
      length = row['Length']
      cat = row['Category']
      if cat == 'ins':
        gapped_aligns.append(trim_alignment(gt, cutsite, 'ins'))
        continue
      if gt_pos == 'e':
        gapped_aligns.append('multiple deletion genotypes')
        continue

      gt_pos = int(gt_pos)
      gap_gt = gt[:cutsite - length + gt_pos] + '-'*length + gt[cutsite - length + gt_pos:]
      gap_gt = add_bar(gap_gt, cutsite)
      gapped_aligns.append(trim_alignment(gap_gt, cutsite, 'del'))
    return gapped_aligns

  gap_gts = get_gapped_alignments(top10, stats)
  
  alignments, categories = [], []
  cutsite = stats['Cutsite'].iloc[0]
  reference_seq = stats['Reference sequence'].iloc[0]
  reference_seq = add_bar(reference_seq, cutsite)
  alignments.append(trim_alignment(reference_seq, cutsite, 'ref'))
  categories.append('Reference')

  for gt, length, cat in zip(gap_gts, lens, cats):
    alignments.append(gt)
    if cat == 'ins':
      categories.append('%s-bp insertion' % (length))
    elif cat == 'del':
      categories.append('%s-bp deletion' % (length))

  return dict(
    data = [go.Table(
      type = 'table',
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

@app.callback(
  Output('S_summary-alignment-barchart', 'figure'),
  [Input('S_hidden-pred-df', 'children'),
   Input('S_hidden-pred-stats', 'children'),
  ])
def update_summary_alignment_barchart(pred_df_string, pred_stats_string):
  pred_df = pd.read_csv(StringIO(pred_df_string), index_col = 0)
  stats = pd.read_csv(StringIO(pred_stats_string), index_col = 0)
  
  inDelphi.add_genotype_column(pred_df, stats)

  top10 = pred_df.sort_values('Predicted frequency', ascending = False).iloc[:10]
  fqs = top10['Predicted frequency'][::-1]

  cats = top10['Category'][::-1]
  gts = top10['Genotype position'][::-1]
  colors = []
  for cat, gt in zip(cats, gts):
    if gt == 'e':
      colors.append('rgb(236, 100, 12)')
    else:
      if cat == 'del':
        colors.append('rgb(221, 46, 31)')
      else:
        colors.append('rgb(0, 160, 220)')

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
  [Input('S_hidden-pred-df', 'children'),
   Input('S_hidden-pred-stats', 'children'),
  ])
def plot_genstats_precision(pred_df_string, pred_stats_string):
  pred_df = pd.read_csv(StringIO(pred_df_string), index_col = 0)
  stats = pd.read_csv(StringIO(pred_stats_string), index_col = 0)
  xval = stats['Precision'].iloc[0]
  return dict(
    data = [
      generalStats.gs_precision.trace(xval),
    ],
    layout = generalStats.gs_precision.layout(xval),
  )

@app.callback(
  Output('S_plot-genstats-logphi', 'figure'),
  [Input('S_hidden-pred-df', 'children'),
   Input('S_hidden-pred-stats', 'children'),
  ])
def plot_genstats_logphi(pred_df_string, pred_stats_string):
  pred_df = pd.read_csv(StringIO(pred_df_string), index_col = 0)
  stats = pd.read_csv(StringIO(pred_stats_string), index_col = 0)
  xval = np.log(stats['Phi'].iloc[0])
  return dict(
    data = [
      generalStats.gs_logphi.trace(xval),
    ],
    layout = generalStats.gs_logphi.layout(xval),
  )

@app.callback(
  Output('S_plot-genstats-frameshift', 'figure'),
  [Input('S_hidden-pred-df', 'children'),
   Input('S_hidden-pred-stats', 'children'),
  ])
def plot_genstats_frameshift(pred_df_string, pred_stats_string):
  pred_df = pd.read_csv(StringIO(pred_df_string), index_col = 0)
  stats = pd.read_csv(StringIO(pred_stats_string), index_col = 0)
  xval = stats['Frameshift frequency'].iloc[0]
  return dict(
    data = [
      generalStats.gs_frameshift.trace(xval),
    ],
    layout = generalStats.gs_frameshift.layout(xval),
  )


## General stats text
@app.callback(
  Output('S_text-genstats-precision', 'children'),
  [Input('S_hidden-pred-df', 'children'),
   Input('S_hidden-pred-stats', 'children'),
  ])
def text_genstats_precision(pred_df_string, pred_stats_string):
  pred_df = pd.read_csv(StringIO(pred_df_string), index_col = 0)
  stats = pd.read_csv(StringIO(pred_stats_string), index_col = 0)
  xval = stats['Precision'].iloc[0]
  cum, var_text, var_color = generalStats.gs_precision.cumulative(xval)
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
  [Input('S_hidden-pred-df', 'children'),
   Input('S_hidden-pred-stats', 'children'),
  ])
def text_genstats_logphi(pred_df_string, pred_stats_string):
  pred_df = pd.read_csv(StringIO(pred_df_string), index_col = 0)
  stats = pd.read_csv(StringIO(pred_stats_string), index_col = 0)
  xval = np.log(stats['Phi'].iloc[0])
  cum, var_text, var_color = generalStats.gs_logphi.cumulative(xval)
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
    html.Span('Log phi: %.2f' % (xval),
      className = 'generalstats_subtext_style'),
    html.Br(),
    html.Span('Percentile: %s' % (cum),
      className = 'generalstats_subtext_style'),
  ]

@app.callback(
  Output('S_text-genstats-frameshift', 'children'),
  [Input('S_hidden-pred-df', 'children'),
   Input('S_hidden-pred-stats', 'children'),
  ])
def text_genstats_frameshift(pred_df_string, pred_stats_string):
  pred_df = pd.read_csv(StringIO(pred_df_string), index_col = 0)
  stats = pd.read_csv(StringIO(pred_stats_string), index_col = 0)
  xval = stats['Frameshift frequency'].iloc[0]
  cum, var_text, var_color = generalStats.gs_frameshift.cumulative(xval)
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
  [Input('S_hidden-pred-df', 'children'),
  ])
def plot_indel_len(pred_df_string):
  pred_df = pd.read_csv(StringIO(pred_df_string), index_col = 0)

  lendf = inDelphi.get_indel_length_fqs(pred_df)

  X = [int(s) for s in lendf['Indel length']]
  Y = [s for s in lendf['Predicted frequency']]
  return dict(
    data = [
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
    ],
    layout = go.Layout(
      xaxis = dict(
        autorange = 'reversed',
        title = 'Indel length',
        ticks = 'outside',
        ticklen = 3,
        tickwidth = 0.5,
        tick0 = 0,
        dtick = 5,
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
  [Input('S_hidden-pred-df', 'children'),
  ])
def plot_fs(pred_df_string):
  pred_df = pd.read_csv(StringIO(pred_df_string), index_col = 0)

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
# Genotype table callbacks
## 
@app.callback(
  Output('S_table-genotypes', 'rows'), 
  [Input('S_hidden-pred-df', 'children'),
   Input('S_hidden-pred-stats', 'children'),
  ])
def update_genotype_table(pred_df_string, pred_stats_string):
  pred_df = pd.read_csv(StringIO(pred_df_string), index_col = 0)
  stats = pd.read_csv(StringIO(pred_stats_string), index_col = 0)

  inDelphi.add_genotype_column(pred_df, stats)
  inDelphi.add_name_column(pred_df, stats)

  # Edit for display
  pred_df['Frequency (%)'] = [float('%.1f' % (s)) for s in pred_df['Predicted frequency']]
  pred_df = pred_df.drop(
    [
      'Inserted Bases', 
      'Genotype position',
      'Predicted frequency',
      'Category',
    ], 
    axis = 1
  )
  pred_df = pred_df[['Name', 'Length', 'Frequency (%)', 'Genotype']]
  return pred_df.to_dict('records')

@app.callback(
  Output('S_plot-table-genotypes', 'figure'),
  [Input('S_table-genotypes', 'rows'),
   Input('S_table-genotypes', 'selected_row_indices')
  ])
def update_genotype_plot(rows, selected_row_indices):
  df = pd.DataFrame(rows)
  colors =  ['#0074D9'] * len(df)
  for idx in (selected_row_indices or []):
    colors[idx] = '#FF851B'
  return dict(
    data = [
      go.Bar(
        x = df['Name'],
        y = df['Frequency (%)'],
        marker = dict(
          color = colors,
          line = dict(
            width = 0,
          ),
        ),
      ),
    ],
    layout = go.Layout(
      xaxis = dict(
        fixedrange = True,
      ),
      yaxis = dict(
        title = 'Frequency (%)',
        fixedrange = True,
      ),
      font = dict(
        family = 'Arial',
      ),
      margin = dict(
        t = 10,
        r = 20,
      ),
    ),
  )

@app.callback(
  Output('S_table-genotypes', 'selected_row_indices'),
  [Input('S_plot-table-genotypes', 'clickData')],
  [State('S_table-genotypes', 'selected_row_indices')]
  )
def update_datatable_selected(clickData, selected_row_indices):
  # Update selections in table based on clicking plot
  if clickData:
    for point in clickData['points']:
      if point['pointNumber'] in selected_row_indices:
        selected_row_indices.remove(point['pointNumber'])
      else:
        selected_row_indices.append(point['pointNumber'])
  return selected_row_indices


##
# Download callbacks
##
@app.callback(
  Output('S_csv-download-link', 'href'), 
  [Input('S_hidden-pred-df', 'children'),
   Input('S_hidden-pred-stats', 'children'),
  ])
def update_link(pred_df_string, pred_stats_string):
  pred_df = pd.read_csv(StringIO(pred_df_string), index_col = 0)
  stats = pd.read_csv(StringIO(pred_stats_string), index_col = 0)

  inDelphi.add_genotype_column(pred_df, stats)
  inDelphi.add_name_column(pred_df, stats)

  time = str(datetime.datetime.now()).replace(' ', '_').replace(':', '-')
  link_fn = '/dash/urlToDownload?value={}'.format(time)
  pred_df.to_csv('user-csvs/%s.csv' % (time))
  return link_fn

@app.callback(
  Output('S_summary-download-link', 'href'),
  [Input('S_hidden-pred-df', 'children'),
   Input('S_hidden-pred-stats', 'children')],
  [State('S_page-link', 'href')])
def update_summary_link(pred_df_string, pred_stats_string, pagelink):
  pred_df = pd.read_csv(StringIO(pred_df_string), index_col = 0)
  stats = pd.read_csv(StringIO(pred_stats_string), index_col = 0)

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
  [Input('S_textbox1', 'value',),
   Input('S_textbox2', 'value',),
  ])
def update_pagelink(text1, text2):
  seq = text1 + text2
  cutsite = len(text1)
  return 'https://dev.crisprindelphi.design%s' % (lib.encode_dna_to_url_path_single(seq, cutsite))

