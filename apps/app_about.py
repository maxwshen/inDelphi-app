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
                dcc.Markdown(dedent('''
                  This website provides an online interactive implementation of inDelphi, a computational method that predicts DNA repair outcomes at DNA double strand breaks induced by CRISPR/SpCas9 resulting from non-homologous end-joining. inDelphi is described in:

                  Max W. Shen*, Mandana Arbab*, Jonathan Y. Hsu, Daniel Worstell, Sannie J. Culbertson, Olga Krabbe, Christopher A. Cassa, David R. Liu, David K. Gifford, and Richard I. Sherwood. "Predictable and precise template-free editing of pathogenic variants." _Nature_, in press, 2018.

                  This user guide page is intended as a resource for biologists using CRISPR. For a layperson-friendly description, please refer to our About page.
                  '''),
                  className = 'markdown_style',
                ),
              ],
              id = 'overview',
              className = 'hashbuffer',
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
