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
  # Hidden divs for light data storage
  ###################################################
  html.Div(
    [
      html.Div(
        id = 'GUIDE_hidden-counter',
        children = 'init'
      ),
      dcc.Location(
        id = 'GUIDE_url',
        refresh = False,
      ),
    ],
    style = dict(
      display = 'none',
    ),
  ),

  ###################################################
  # Header
  ###################################################
  html.Div(
    [
      ###################################################
      # Upper header
      ###################################################
      header.get_navigation_header('guide'),

      html.Div(
        style = dict(
          borderBottom = '3px solid #777777',
          width = '1010px',
          margin = '0 auto',
        ),
      )

    ],
    # style = dict(
    #   position = 'fixed',
    #   top = 0,
    #   backgroundColor = 'white',
    #   width = '1010px',
    #   left = '50%',
    #   zIndex = 10,
    #   transform = 'translate(-50%, 0px)',
    #   borderBottom = '3px solid #777777',
    #   marginBottom = '30px',
    # ),
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

  # Side bar
  html.Div(
    [


      ## 
      html.Div([
        html.A('Overview', 
          href = 'guide#overview', 
          className = 'dynamicunderline',
          style = dict(color = 'black', textDecoration = 'none')),
      ], style = dict(marginBottom = '8px')),

      html.Div([
        html.A('Description', 
          href = 'guide#overview2', 
          className = 'dynamicunderline',
          style = dict(color = 'black', textDecoration = 'none')),
      ], style = dict(paddingLeft = '24px', marginBottom = '8px')),

      html.Div([
        html.A('Overview of models', 
          href = 'guide#overview3', 
          className = 'dynamicunderline',
          style = dict(color = 'black', textDecoration = 'none')),
      ], style = dict(paddingLeft = '24px', marginBottom = '8px')),

      html.Div([
        html.A('General website features', 
          href = 'guide#overview4', 
          className = 'dynamicunderline',
          style = dict(color = 'black', textDecoration = 'none')),
      ], style = dict(paddingLeft = '24px', marginBottom = '8px')),

      html.Div([
        html.A('Python implementation', 
          href = 'guide#overview5', 
          className = 'dynamicunderline',
          style = dict(color = 'black', textDecoration = 'none')),
      ], style = dict(paddingLeft = '24px', marginBottom = '8px')),

      ##
      html.Div([
        html.A('Single mode & inDelphi model details', 
          href = 'guide#single', 
          className = 'dynamicunderline',
          style = dict(color = 'black', textDecoration = 'none')),
      ], style = dict(marginBottom = '8px')),

      html.Div([
        html.A('inDelphi predicts 80-95% of editing outcomes', 
          href = 'guide#single', 
          className = 'dynamicunderline',
          style = dict(color = 'black', textDecoration = 'none')),
      ], style = dict(paddingLeft = '24px', marginBottom = '8px')),

      html.Div([
        html.A('Microhomology and overlap length inform deletions', 
          href = 'guide#single2', 
          className = 'dynamicunderline',
          style = dict(color = 'black', textDecoration = 'none')),
      ], style = dict(paddingLeft = '24px', marginBottom = '8px')),

      html.Div([
        html.A('1-bp insertions depend on cell-type, DNA motif, and gRNA orientation', 
          href = 'guide#single3', 
          className = 'dynamicunderline',
          style = dict(color = 'black', textDecoration = 'none')),
      ], style = dict(paddingLeft = '24px', marginBottom = '8px')),

      html.Div([
        html.A('Microhomology-less deletions are noisy', 
          href = 'guide#single4', 
          className = 'dynamicunderline',
          style = dict(color = 'black', textDecoration = 'none')),
      ], style = dict(paddingLeft = '24px', marginBottom = '8px')),

      html.Div([
        html.A('Features: Single mode', 
          href = 'guide#single5', 
          className = 'dynamicunderline',
          style = dict(color = 'black', textDecoration = 'none')),
      ], style = dict(paddingLeft = '24px', marginBottom = '8px')),

      ##
      html.Div([
        html.A('Batch mode', 
          href = 'guide#batch', 
          className = 'dynamicunderline',
          style = dict(color = 'black', textDecoration = 'none')),
      ], style = dict(marginBottom = '8px')),

      html.Div([
        html.A('Features', 
          href = 'guide#batch2', 
          className = 'dynamicunderline',
          style = dict(color = 'black', textDecoration = 'none')),
      ], style = dict(paddingLeft = '24px', marginBottom = '8px')),

      html.Div([
        html.A('Processing limitations', 
          href = 'guide#batch3', 
          className = 'dynamicunderline',
          style = dict(color = 'black', textDecoration = 'none')),
      ], style = dict(paddingLeft = '24px', marginBottom = '8px')),

      html.Div([
        html.A('Summary statistics', 
          href = 'guide#batch4', 
          className = 'dynamicunderline',
          style = dict(color = 'black', textDecoration = 'none')),
      ], style = dict(paddingLeft = '24px', marginBottom = '8px')),

      ##
      html.Div([
        html.A('Gene mode', 
          href = 'guide#gene', 
          className = 'dynamicunderline',
          style = dict(color = 'black', textDecoration = 'none')),
      ], style = dict(marginBottom = '8px')),

      html.Div([
        html.A('Features', 
          href = 'guide#gene2', 
          className = 'dynamicunderline',
          style = dict(color = 'black', textDecoration = 'none')),
      ], style = dict(paddingLeft = '24px', marginBottom = '8px')),

      html.Div([
        html.A('Summary statistics', 
          href = 'guide#gene3', 
          className = 'dynamicunderline',
          style = dict(color = 'black', textDecoration = 'none')),
      ], style = dict(paddingLeft = '24px', marginBottom = '8px')),

      ##
      html.Div([
        html.A('FAQ', 
          href = 'guide#faq', 
          className = 'dynamicunderline',
          style = dict(color = 'black', textDecoration = 'none')),
      ], style = dict(marginBottom = '8px')),


    ],
    style = dict(
      backgroundColor = 'white',
      opacity = 0.9,
      zIndex = 20,
      top = '15%',
      left = '2%',
      overflowX = 'hidden',
      position = 'fixed',
      width = '15%',
      # height = 300,
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
                  ## Overview
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

            html.Div(
              [
                dcc.Markdown(dedent('''
                  #### Description
                  inDelphi is a predictive model for events following CRISPR cleavage at a specified target site that result in a non-wild type genotype detectable by Illumina sequencing. inDelphi predicts the relative frequencies of 1-bp insertion and 1- to 60-bp deletion genotypes which comprise 80-95% of editing outcomes observed in our data. inDelphi can be used, for instance, to predict the frameshift frequency of a gRNA targeting a coding region of a gene. 

                  inDelphi is unconcerned with off-target editing, where a single gRNA may induce cleavage at multiple target sites. inDelphi is also unconcerned with the likelihood of a cleavage event, a.k.a. cutting efficiency (addressed, for example, by Azimuth: [Doench, Fusi et al., 2016](https://www.nature.com/articles/nbt.3437)), or frequencies of repair back to wild-type. inDelphi is trained on data with long exposure to CRISPR with saturated levels of indels, and likely is not applicable to experiments with brief exposure to CRISPR. inDelphi was not trained on data sensitive to large (1kb+) deletions.

                  The inDelphi model was trained on data in mESC, U2OS, HCT116, K562, and HEK293 cell types generated by ourselves and [van Overbeek et al., 2016](https://www.cell.com/molecular-cell/fulltext/S1097-2765(16)30325-2). While we are most confident of predictions in these cell-types, and though we observe moderate to large cell-type variability in our data, we do expect inDelphi to be relevant beyond these five cell types to other mammalian cell types. Human embryonic stem cells, for instance, are likely to have similar repair outcomes as mESCs. inDelphi is not expected to generalize well to bacteria, plants, and non-mammalian eukaryotes such as yeast.

                  inDelphi was trained for SpCas9. With KKH-SaCas9, we currently hypothesize that predictions of deletions are relevant while 1-bp insertions are not. inDelphi is not expected to work for Cpf1/Cas12a.

                  inDelphi was trained and validated on both artificial genome-integrated and natively endogenous DNA target sites, ranging in GC content, microhohomology strength, heterogeneity of indel outcomes, and Azimuth-predicted cutting efficiency. We expect it to be relevant to native endogenous DNA target sites.

                  '''),
                  className = 'markdown_style',
                ),
              ],
              id = 'overview2',
              className = 'hashbuffer',
            ),


            html.Div(
              [
                dcc.Markdown(dedent('''
                  #### Overview of modes

                  This online interactive implementation of inDelphi supports:

                  **1. Single mode**

                  Input: A sequence context and cleavage location.

                  Output: Indel frequency predictions and summary statistics.
                  
                  **2. Batch mode**

                  Input: A sequence context and PAM.

                  Output: Summary statistics of predictions at all gRNAs in the sequence context compatible with the PAM. Our online tool is limited to 80 gRNAs.
                  
                  **3. Gene mode**

                  Input: A hg38 or mm10 coding gene.

                  Output: Summary statistics of predictions at all SpCas9 gRNAs targeting the coding gene.
                  '''),
                  className = 'markdown_style',
                ),
              ],
              id = 'overview3',
              className = 'hashbuffer',
            ),

            html.Div(
              [
                dcc.Markdown(dedent('''
                  #### General website features
                  - Interactivity! Created with Dash and plotly
                  - Five cell-type specific versions of inDelphi
                  - Single mode, batch mode, gene mode
                  - Download tables of all predictions
                  - Shareable URLs save your settings
                  '''),
                  className = 'markdown_style',
                ),
              ],
              id = 'overview4',
              className = 'hashbuffer',
            ),

            html.Div(
              [
                dcc.Markdown(dedent('''
                  #### Python implementation

                  A Python implementation of the inDelphi model is available in this [github repository](https://github.com/maxwshen/inDelphi-model), and is recommended for computational users interested in running inDelphi at a larger scope than supported by our online implementation.

                  '''),
                  className = 'markdown_style',
                ),
              ],
              id = 'overview5',
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


    ###################################################
    # Section: Single
    ###################################################
    html.Div(
      [
        html.Div(
          [
            dcc.Markdown(dedent('''
              ## Single mode & inDelphi model details 

              #### inDelphi predicts 80-95% of editing outcomes
              We interrogated indel repair at thousands of designed genome-integrated target sites and ~90 endogenous DNA target sites with SpCas9 in mESCs, U2OS, HEK293, HCT116, and K562 cells. Each target site yielded dozens to hundreds of unique indels. Across cell-types, we observe that three repair classes constitute 80-95% of repair outcomes: 1-bp insertions, microhomology deletions, and microhomology-less deletions (Fig. 1). inDelphi models these three repair classes.
              '''),
              className = 'markdown_style',
            ),

            html.Div(
              [
                html.Img(src = '/assets/fig-pies.PNG'),
                html.Div(
                  [
                    'Fig. 1: Observed breakdown of DNA repair types by cell-type. MH: microhomology. Endo.: endogenous.'
                    # 'Figure 2: inDelphi predicts +1 bp insertions and 1- to 60-bp deletions (negative indel length) separated into microhomology deletions (red) and microhomology-less deletions (orange).'
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
              Microhomology deletions are defined as deletions that are consistent with sequence microhomology. inDelphi recognizes and includes in its modeling the fact that, under this definition, microhomology deletions can arise from both microhomology-dependent repair mechanisms like microhomology-mediated end-joining (MMEJ) and microhomology-independent mechanisms like Lig4-mediated directed end-joining.

              The remaining repair outcomes are rare and comprise 2+ bp insertions, indels occurring nearby but not at the target site which are enriched in Cas9 treatment vs. control cells, and combination indels. inDelphi does not model these classes of repair due to being relatively rare and having higher noise levels (lower consistency between experimental replicates).

              inDelphi is trained on Illumina deep sequencing data that is unable to detect long (1kb+) deletions. inDelphi models 1- to 60-bp deletions, which comprise the vast majority of deletions observed with Illumina sequencing.
              '''),
              className = 'markdown_style',
            ),

            ###################################################
            # Subsection
            ###################################################
            html.Div(
              [
                dcc.Markdown(dedent('''
                  #### Microhomology and overlap length inform deletions
                  '''),
                  className = 'markdown_style',
                ),
              ],
              id = 'single2',
              className = 'hashbuffer',
            ),

            html.Div(
              [
                html.Img(src = '/assets/fig-indel-len.gif'),
                html.Div(
                  [
                    'Fig. 2: Example inDelphi predictions. inDelphi predicts +1 bp insertions and 1- to 60-bp deletions (negative indel length) separated into microhomology deletions (red) and microhomology-less deletions (orange).'
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
              Microhomology-less deletions (orange) decay in frequency with increasing deletion length. Microhomology deletions also decay with increasing deletion length, but strong sequence microhomology can compensate. We model microhomology deletions with the microhomology-mediated end-joining mechanism (Fig. 3). 
              '''),
              className = 'markdown_style',
            ),

            html.Div(
              [
                html.Img(src = '/assets/fig-mechanism.PNG'),
                html.Div(
                  [
                    'Fig. 3: DNA repair mechanisms leading to microhomology and microhomology-less deletions. Left: microhomology-mediated end-joining (MMEJ). Right: canonical non-homologous end-joining (c-NHEJ) mediated by Lig4.'
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

            ###################################################
            # Subsection
            ###################################################
            html.Div(
              [
                dcc.Markdown(dedent('''
                  #### 1-bp insertions depend on cell-type, DNA motif, and gRNA orientation
                  '''),
                  className = 'markdown_style',
                ),
              ],
              id = 'single3',
              className = 'hashbuffer',
            ),

            html.Div(
              [
                html.Img(src = '/assets/fig-1bpins.PNG'),
                html.Div(
                  [
                    'Fig. 4: 1-bp insertions copy the 5\' adjacent nucleotide (a), and occur at a frequency related to a DNA motif in mESCs (b) and U2OS cells (c).'
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
              We observe that the major axis of cell-type variability is in the relative frequency of 1-bp insertions to deletions (Fig. 1, reproduced below. Compare the size of the blue arcs across cell-types). 1-bp insertions copy the nucleotide at position -4, which we label as the PAM-distal nucleotide that is directly adjacent to the cleavage site (Fig. 4a). Despite this cell-type variability, a single DNA motif emerges in mESC and U2OS cells when learning a predictive model of 1-bp insertion frequency (Fig. 4b,c).

              inDelphi uses the same computational model for deletions across all cell-types, only using cell-type information when predicting the relative frequency of 1-bp insertions to deletions.

              1-bp insertions appear to be an SpCas9-specific phenomena. Experiments by [Lemos et al., 2018](http://www.pnas.org/content/early/2018/02/12/1716855115) demonstrate that two SpCas9 gRNAs inducing the same cleavage site from both sequence orientations consistently copy the PAM-distal nucleotide. Preliminary unpublished experiments by us indicate that KKH-SaCas9 does not share the same 1-bp insertion pattern.
              '''),
              className = 'markdown_style',
            ),

            html.Div(
              [
                html.Img(src = '/assets/fig-pies.PNG'),
                html.Div(
                  [
                    'Fig. 1: Observed breakdown of DNA repair types by cell-type. MH: microhomology. Endo.: endogenous.'
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

            ###################################################
            # Subsection
            ###################################################
            html.Div(
              [
                dcc.Markdown(dedent('''
                  #### Microhomology-less deletions are noisy
                  '''),
                  className = 'markdown_style',
                ),
              ],
              id = 'single4',
              className = 'hashbuffer',
            ),

            html.Div(
              [
                html.Img(src = '/assets/fig-mhless.png'),
                html.Div(
                  [
                    'Fig. 5: Example of separating microhomology-less deletions of a single deletion length into genotype-resolution predictions.'
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
              inDelphi, as presented in our manuscript, is not evaluated on its performance of predicting microhomology-less deletions at single-base resolution. In our online implementation, we provide single-base resolution predictions for user convenience, and highlight caveats and warnings here.

              Microhomology-less deletions are the most noisy of the three repair classes modeled by inDelphi, with the lowest consistency between experimental replicates. They constitute about 30-40% of all deletions, and their most obvious pattern is an exponentially decreasing frequency with increasing deletion length.

              In our manuscript, we evaluate predictive performance on microhomology-less deletions at deletion length resolution. In the example in Fig. 5, there are 4 unique genotypes corresponding to a microhomology-less deletion of length 3. inDelphi predicts the sum total frequency of all 4 genotypes together: 4.86%. 

              To achieve single-base resolution predictions, we perform an additional unvalidated heuristic of attributing 1/3 to each unilateral end-joining possibility (total 2/3, reflecting the dominance of Lig4-mediated end-joining, refer to the DNA repair mechanisms in Fig. 3), and the last 1/3 is split evenly across all other possibilities.

              We do not recommend relying on microhomology-less deletions for single-base resolution applications, such as correcting pathogenic alleles to wild type genotypes. We recommend putting more faith into microhomology deletions and SpCas9 1-bp insertions.

              '''),
              className = 'markdown_style',
            ),


            ###################################################
            # Subsection
            ###################################################
            html.Div(
              [
                dcc.Markdown(dedent('''
                  #### Features: Single mode

                  Input: A sequence context and cleavage location.

                  Output: Indel frequency predictions and summary statistics.

                  - IUPAC PAMs are supported, such as NGG, NNNRRT
                  - Calculates summary statistics such as frameshift frequency, indel length frequency, precision, microhomology strength
                  - Comparison of summary statistics of your result to precomputed statistics at 13 million SpCas9 gRNAs targeting human exons and introns
                  - Shareable URL that saves your specified target site and cut site location
                  '''),
                  className = 'markdown_style',
                ),
              ],
              id = 'single5',
              className = 'hashbuffer',
            ),

          ],
          style = dict(
            width = '800px',
            margin = '0 auto',
          ),
        ),


      ],
      id = 'single',
      className = 'hashbuffer',
    ),

    ###################################################
    # Section: Batch
    ###################################################
    html.Div(
      [
        html.Div(
          [
            dcc.Markdown(dedent('''
              ## Batch mode

              '''),
              className = 'markdown_style',
            ),

            html.Div(
              [
                html.Img(src = '/assets/fig-batch.gif'),
                html.Div(
                  [
                    'Fig. 6: Example usage of batch mode, showcasing gRNA selection and sorting.'
                  ],
                  style = dict(
                    fontStyle = 'italic',
                    width = '450px',
                    margin = '0 auto',
                    marginBottom = '20px',
                    transform = 'translateX(212px)',
                  )
                ),
              ],
              style = dict(
                textAlign = 'center',
                transform = 'translateX(-212px)',
              ),
            ),

            ###################################################
            # Subsection
            ###################################################
            html.Div(
              [
                dcc.Markdown(dedent('''
                  #### Features

                  Input: A sequence context and PAM.

                  Output: Summary statistics of predictions at all gRNAs in the sequence context compatible with the PAM. Our online tool is limited to 80 gRNAs.

                  Note: This tool will only consider gRNAs for which at least 30-bp of sequence exists on both sides of the cleavage site. This is because inDelphi requires sequence context surrounding the cleavage site to make predictions.

                  - IUPAC PAMs are supported, such as NGG, NNNRRT
                  - Can score gRNAs for their frequency of repairing to a specific genotype (e.g., correcting a pathogenic allele), or deleting particular nucleotides (e.g., ablating a splite site), or their distance to a position of interest (e.g., selecting gRNAs with high frameshift frequency that are deep inside an exon and far from an exon boundary)
                  - Calculates many summary statistics (more details below)
                  - Links for all gRNAs to single mode, enabling genotype-resolution inspection 
                  - Shareable URL that saves your specified target site, PAM, advanced options, displayed columns, sort by column and direction, and selected gRNA
                  - Select gRNAs by clicking for visualization and comparison
                  '''),
                  className = 'markdown_style',
                ),
              ],
              id = 'batch2',
              className = 'hashbuffer',
            ),

            ###################################################
            # Subsection
            ###################################################
            html.Div(
              [
                dcc.Markdown(dedent('''
                  #### Processing limitations

                  The online implementation is limited to 80 gRNAs. This is to avoid requests taking longer than 30 seconds, which our hosting service Heroku will time out and fail. If you want to run inDelphi at a larger scope, you can download our Python implementation and run it locally or on your server.
                  '''),
                  className = 'markdown_style',
                ),
              ],
              id = 'batch3',
              className = 'hashbuffer',
            ),

            ###################################################
            # Subsection
            ###################################################
            html.Div(
              [
                dcc.Markdown(dedent('''
                  #### Summary statistics

                  Cutsite: 0-based index of the cutsite within your provided sequence.

                  Precision: A statistical measure (negative entropy normalized to range between 0 and 1) of the frequency distribution of indels for a gRNA. High precision is closer to 1 and represents gRNAs where a small number of repair outcomes constitute nearly all outcomes. Low precision is closer to 0 and represents gRNAs with no dominant repair outcomes.

                  Frameshift (%): Total predicted frequency of frame +1 and +2 indels out of all 1-bp insertions and 1- to 60-bp deletions, which in total constitute 80-95% of repair outcomes. Does not consider wild-type outcomes. The typical frameshift % is higher than 2/3 because 1-bp insertions and 1- to 2-bp deletions are enriched.

                  Frame +0, +1, +2 (%): Total predicted frequency of indels resulting in that frame change out of all 1-bp insertions and 1- to 60-bp deletions.

                  Microhomology strength (MH strength): A relative score representing the sum total strength of all sequence microhomologies. Longer microhomology length and higher GC content improves strength, while deletion length decreases strength. Scale is arbitrary: 0 has no particular meaning, neither does negative values.

                  Most frequent genotype / deletion / insertion (%) (M.F. gt / del / ins (%)): The predicted frequency of the most frequent genotype. Related to precision which considers the full frequency distribution, but this statistic only considers the single most frequent genotype.

                  Expected indel length (Exp. indel len): The average indel length out of all 1-bp insertions and 1- to 60-bp deletions. Does not consider wild-type outcomes. Positive numbers here indicate deletion length.

                  '''),
                  className = 'markdown_style',
                ),
              ],
              id = 'batch4',
              className = 'hashbuffer',
            ),

          ],
          style = dict(
            width = '800px',
            margin = '0 auto',
          ),
        ),
      ],
      id = 'batch',
      className = 'hashbuffer',
    ),


    ###################################################
    # Section: Gene
    ###################################################
    html.Div(
      [
        html.Div(
          [
            dcc.Markdown(dedent('''
              ## Gene mode

              '''),
              className = 'markdown_style',
            ),

            html.Div(
              [
                html.Img(src = '/assets/fig-gene.gif'),
                html.Div(
                  [
                    'Fig. 7: Example usage of gene mode on IL32 in hg38.'
                  ],
                  style = dict(
                    fontStyle = 'italic',
                    width = '450px',
                    margin = '0 auto',
                    marginBottom = '20px',
                    transform = 'translateX(300px)',
                  )
                ),
              ],
              style = dict(
                textAlign = 'center',
                transform = 'translateX(-300px)',
              ),
            ),

            ###################################################
            # Subsection
            ###################################################
            html.Div(
              [
                dcc.Markdown(dedent('''
                  #### Features

                  Input: A hg38 or mm10 coding gene.

                  Output: Summary statistics of predictions at all SpCas9 gRNAs targeting the coding gene.

                  Note: The online visualization is limited to viewing kGIDs that have a sum total of fewer than 1000 gRNAs. You can, however, download the full table of predictions.

                  - Responsive search for coding genes. Note: Query is matched anywhere in gene name, not necessarily from the beginning. E.g., "Mt1" will match Dnmt1. You can scroll to find your desired gene.
                  - Supports all kgIDs described in the UCSC genome browser
                  - Calculates many summary statistics (more details below), including exon number, and distance to 5' and 3' exon boundaries
                  '''),
                  className = 'markdown_style',
                ),
              ],
              id = 'gene2',
              className = 'hashbuffer',
            ),

            ###################################################
            # Subsection
            ###################################################
            html.Div(
              [
                dcc.Markdown(dedent('''
                  #### Summary statistics

                  Exon number: The exon number for that kgID.

                  Distance to 5'/3' exon boundary (Dist. to 5'/3' end): Distance of the cutsite to the 5' or 3' exon boundary. gRNAs causing indels near an exon boundary may influence splicing or cause unintended consequences instead of, or beyond, knocking out a coding gene.

                  Precision: A statistical measure (negative entropy normalized to range between 0 and 1) of the frequency distribution of indels for a gRNA. High precision is closer to 1 and represents gRNAs where a small number of repair outcomes constitute nearly all outcomes. Low precision is closer to 0 and represents gRNAs with no dominant repair outcomes.

                  Frameshift (%): Total predicted frequency of frame +1 and +2 indels out of all 1-bp insertions and 1- to 60-bp deletions, which in total constitute 80-95% of repair outcomes. Does not consider wild-type outcomes. The typical frameshift % is higher than 2/3 because 1-bp insertions and 1- to 2-bp deletions are enriched.

                  Frame +0, +1, +2 (%): Total predicted frequency of indels resulting in that frame change out of all 1-bp insertions and 1- to 60-bp deletions.

                  Microhomology strength (MH strength): A relative score representing the sum total strength of all sequence microhomologies. Longer microhomology length and higher GC content improves strength, while deletion length decreases strength. Scale is arbitrary: 0 has no particular meaning, neither does negative values.

                  Most frequent genotype / deletion / insertion (%) (M.F. gt / del / ins (%)): The predicted frequency of the most frequent genotype. Related to precision which considers the full frequency distribution, but this statistic only considers the single most frequent genotype.

                  Expected indel length (Exp. indel len): The average indel length out of all 1-bp insertions and 1- to 60-bp deletions. Does not consider wild-type outcomes. Positive numbers here indicate deletion length.

                  '''),
                  className = 'markdown_style',
                ),
              ],
              id = 'gene3',
              className = 'hashbuffer',
            ),

          ],
          style = dict(
            width = '800px',
            margin = '0 auto',
          ),
        ),
      ],
      id = 'gene',
      className = 'hashbuffer',
    ),

    ###################################################
    # Section: FAQ
    ###################################################
    html.Div(
      [
        html.Div(
          [
            dcc.Markdown(dedent('''
              ## FAQ and general help

              #### How do I pick a cell-type version of inDelphi?
              We recommend using mESC if your cell-type does not have major defects in DNA repair. 

              #### The website is slow and/or jobs are not being completed.
              Try again at a time with fewer concurrent users. If a job takes more than 30 seconds, it will automatically fail. Parameters are restricted such that most jobs in most situations should not exceed 30 seconds, but this can occur if there are a large number of concurrent users.

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
      id = 'faq',
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
