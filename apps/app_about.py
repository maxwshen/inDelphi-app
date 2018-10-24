import pickle, copy, os, datetime, subprocess
from collections import defaultdict
import numpy as np
import pandas as pd
from scipy.stats import entropy
import time
from io import StringIO
from textwrap import dedent

import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_table_experiments as dt
from dash.dependencies import Input, Output, State
import plotly.graph_objs as go
import flask
import plotly

import inDelphi
import generalStats
import lib, header

from indelphi_app import app


# Remove these plotly modebar buttons to limit interactivity
modebarbuttons_2d = ['zoom2d', 'pan2d', 'select2d', 'lasso2d', 'zoomIn2d', 'zoomOut2d', 'autoScale2d', 'resetScale2d', 'hoverClosestCartesian', 'hoverCompareCartesian', 'toggleSpikelines']


###################################################################
###################################################################
##
# App layout
##
layout = html.Div([

  ###################################################
  # Header
  ###################################################
  html.Div(
    [
      ###################################################
      # Upper header
      ###################################################
      header.get_navigation_header('about'),

      html.Div(
        style = dict(
          borderBottom = '3px solid #777777',
          width = '1010px',
          margin = '0 auto',
        ),
      )

    ],
    style = dict(
      position = 'fixed',
      top = 0,
      backgroundColor = 'white',
      width = '1500px',
      left = '50%',
      zIndex = 10,
      transform = 'translate(-50%, 0px)',
      marginBottom = '30px',
    ),
  ),

  ###################################################
  # Text
  ###################################################
  html.Div(
    [
    ###################################################
    # Section: Overview
    ###################################################
    html.Div(
      [
        html.Div(
          [
            html.Div(
              [
                html.Img(src = '/assets/fig-indel-len.gif'),
                html.Div(
                  [
                    ''
                  ],
                  style = dict(
                    fontStyle = 'italic',
                    width = '450px',
                    margin = '0 auto',
                    marginBottom = '20px',
                  )
                ),
              ],
              style = dict(
                textAlign = 'center',
              ),
            ),

            dcc.Markdown(dedent('''
              inDelphi is a machine learning algorithm that is aimed to assist scientists using CRISPR. Here, we provide a layperson-friendly description of inDelphi. For a more scientifically accurate description, please refer to the User Guide page and our manuscript:

              Max W. Shen*, Mandana Arbab*, Jonathan Y. Hsu, Daniel Worstell, Sannie J. Culbertson, Olga Krabbe, Christopher A. Cassa, David R. Liu, David K. Gifford, and Richard I. Sherwood. "Predictable and precise template-free editing of pathogenic variants." _Nature_, in press, 2018. 

              #### Briefly...
              '''),
              className = 'markdown_style',
            ),

            html.Div(
              [
                html.Img(src = '/assets/fig-coverplus.PNG'),
                html.Div(
                  [
                    ''
                  ],
                  style = dict(
                    fontStyle = 'italic',
                    width = '450px',
                    margin = '0 auto',
                    marginBottom = '20px',
                  )
                ),
              ],
              style = dict(
                textAlign = 'center',
              ),
            ),

            dcc.Markdown(dedent('''
              CRISPR is a genome editing tool that has achieved widespread use by making it easier and cheaper than ever before to induce a DNA double-strand break at a precise location in the genome, specified by a user-designed guide RNA (gRNA). Double-strand breaks in the genome are traumatic to the cell as it will die if it cannot repair it. Cells repair this DNA double-strand break using a variety of _DNA repair pathways_.

              _Homology-directed repair (HDR)_ is one DNA repair pathway that has been in the spotlight. Following a harmful DNA double-strand break, cells use this pathway to repair the affected DNA sequence by copying from the sister-chromatid, by finding areas of similarity surrounding the break-site and filing in the blanks at the break site. Scientists can use this pathway by introducing a DNA template that looks similar to the DNA surrounding the CRISPR cut while carrying a modification at the break site. The cell will use this DNA template to repair the cut and introduce the modification at the cut. The result is controllable editing of DNA at a precise location in the genome into user-specified DNA.

              It turns out, however, that HDR is often less efficient than another DNA repair pathway called _end-joining_. In contrast to HDR which can be used to introduce user-specified edits, end-joining has not been widely recognized as controllable. In the context of not only a traumatic DNA double-strand break, but CRISPR which will repeatedly cut the DNA, end joining can be thought of as a "willy-nilly", no-holds-barred effort to repair the cut into a form that CRISPR can no longer recut. The results are "ugly" insertions and deletions which are hugely diverse and heterogeneous in different cells. While HDR has been a golden child of genome editing, end-joining has been viewed as undesirable noise that is unfortunately more efficient than HDR.

              inDelphi is a computational model that predicts the heterogeneous (100+ unique) mixture of indels resulting from microhomology-mediated end-joining (MMEJ) and non-homologous end-joining (NHEJ) at a CRISPR-induced cut. inDelphi synthesizes known biological mechanisms of DNA repair with machine learning and achieves strong accuracy.

              In our study, we used inDelphi to identify a previously underappreciated class of disease-related mutations that could be corrected to its common, healthy DNA sequence in the majority of CRISPR repairs using the "noisy" end joining pathways. We believe that our study suggests that, with inDelphi, end joining repair can be reinterpreted from an undesirable yet efficient DNA repair pathway, into an efficient DNA repair pathway (1) whose heterogeneity can be predicted and thereby controlled, and (2) that is useful for certain genome editing applications.

              #### Some applications:
              - Identify gRNAs that have extremely high frameshift frequencies to more efficiently knock down genes
              - Identify gRNAs that precisely and efficiently delete specific sequences, such as splicing elements
              - Identify gRNAs with highly precise indel outcomes, to consolidate the risk of undesired NHEJ editing in an HDR application
              - Identify gRNAs with imprecise indel outcomes to boost the information content of CRISPR scars for lineage tracing
              - Depending on the nearby sequence context, identify gRNAs that can introduce a specific desired indel more efficiently than HDR

              #### Credit
              The inDelphi website was built by Max W. Shen using Python, Dash and Plotly. You can contact him at maxwshen (at) mit.edu.

              '''),
              className = 'markdown_style',
            ),

          ],
          style = dict(
            width = '800px',
            margin = '0 auto',
          ),
        ),
      ],
      id = 'overview',
      className = 'hashbuffer',
    ),



    ],
    style = dict(
      transform = 'translateY(120px)',
      marginBottom = '150px',
    )
  )
  ],  # body div
  style = dict(
    width = '1150px',
    margin = '0 auto',
  )
)

#######################################################################
#########################      CALLBACKS      #########################
#######################################################################

@app.server.route('/assets/<resource>')
def serve_image_assets_about(resource):
  # BE VERY CAREFUL NOT TO SERVE ARBITRARY FILES
  return flask.send_from_directory(os.getcwd() + '/assets/', resource)
