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

            dcc.Markdown(dedent('''
              Limited Copyright License for Research Use by Non-Profit and Government Institutions

              BY DOWNLOADING THE CODE OR USING THE SERVICE AND/OR SOFTWARE APPLICATION ACCOMPANYING THIS LICENSE, YOU ARE CONSENTING TO BE BOUND BY ALL OF THE TERMS OF THIS LICENSE

              “Copyright 2018. Massachusetts Institute of Technology, The Broad Institute, Harvard University and Brigham and Women's Hospital. All Rights Reserved.”

              The software is being provided as a service for research, educational, instructional and non-commercial purposes only. By generating a user account and/or submitting jobs to InDelphi you agree to the terms and conditions herein.

              You are an actively enrolled student, post-doctoral researcher, or faculty member at a degree-granting educational institution or US government research institution; and You will only use the InDelphi Software Application and/or Service for educational, instructional, and/or non-commercial research purposes; 

              You understand that all results produced using the Code may only be used for non-commercial research and/or academic purposes;

              You understand that to obtain any right to use the Code for commercial purposes, or in the context of industrially sponsored research, You must enter into an appropriate, separate and direct license agreement with the Owners.

              You will not redistribute unmodified versions of the Code;

              You will redistribute modifications, if any, under the same terms as this license and only to non-profits and US government institutions;

              You must credit the authors of the Code: David K. Gifford, Jonathan Yee-Ting Hsu and Max Walt Shen and cite "Predictable and precise template-free editing of pathogenic mutations by CRISPR-Cas9 nuclease", Nature, 2018, doi:10.1038/s41586-018-0686-x; and

              You understand that neither the names of the Owners nor the names of the authors may be used to endorse or promote products derived from this software without specific prior written permission.
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