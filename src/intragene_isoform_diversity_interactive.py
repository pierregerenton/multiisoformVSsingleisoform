import pandas as pd
import argparse
from dash import Dash, html, dash_table, dcc, Input, Output, callback


parser = argparse.ArgumentParser(
                    prog='visualize',
                    description='visualize',
                    epilog='For more information, contact pierre.gerenton@crg.eu')


parser.add_argument(
    '-d', '--data',
    type=str, required=True,
    help='Path of the output of the "intragene_isoform_diversity" script'
)


args = parser.parse_args()


MARKDOWN_PRESENTATION = """
## Metrics computed in all isoform :

Let $G$, the set, a set of isoform defined as :
                             
$$G = \{ I_1, I_2, \ldots, I_n \}$$
                             
where each $I_i$ is a set of multiple GO terms defined as :
                             
$$ I_i = \{ T_{i1}, T_{i2}, \ldots, T_{im_i} \} $$


### *Number of isoform : number of isoform*
                             
$$n_{isoform} = n = |G|$$

### *Standard deviation of the number of GO term*
                             
$$\sigma_{m_i} = \sqrt{\\frac{1}{n}\sum_{i=1}^{n}{(m_i-\\bar{m_i} )^2}}$$

### *Redudancy metric*
                             
This metrics was designed to have an idea of the number of times a GO term appear in the genes.
For each unique GO term, his number of reoccurence is count ($0$ if it appear only in $1$ isoform, and $n-1$ if it appear in all isoform). Then, the mean of this counting is done to have the mean count of reoccurence. Finally, this count is divided by the $n-1$.
If $r$ is close to $1$, that's means that all GO terms are present in all isoforms, and if $r$ is close to $0$, each isoform is different.

If there is 1 isoform, the redudancy metric is set to $1$.

Let $O$ be the set of all unique GO terms defined as :
                             
$$O = \\bigcup_{i=1}^n I_i = \{ T_{1}, T_{2}, \ldots, T_{n_o} \}$$
                             
Let $count(T_i)$ the number of isoform where $T_i$ is present.
                             
$$rdd(T_i) = \\frac{1}{n-1} \sum_{i=1}^ {n_o}(count(T_i)-1)$$

We then return the mean $\overline{rdd}$.                             


## Metrics computed in all isoform :
Certain metrics were calculated for each pair of isoforms before all the values were averaged.

### *Jaccard index*
The Jaccard index is measure of similarity.
                             
If $I_1 \cup I_2  = \emptyset$, $J(I_1, I_2) = 1$, else

$$J(I_1, I_2) = \\frac{|I_1 \cap I_2|}{|I_1 \cup I_2|}$$

### *Dice coefficient*
                             
The Sørensen–Dice coefficient is measure of similarity.
                             
If $I_1 \cup I_2  = \emptyset$, $D(I_1, I_2) = 1$, else

$$D(I_1, I_2) = \\frac{2|I_1 \cap I_2|}{|I_1| + |I_2|}$$

### *Overlap coefficient*
                             
If $min(|I_1|,|I_2|) = 0$, then $overlap(I_1, I_2) = 1$, else

$$overlap(I_1, I_2) = \\frac{|I_1 \cap I_2|}{min(|I_1|,|I_2|)}$$

### *BP, CC and MF GOGO similarity*
Semantic similarity computed between GO terms of each isoform.
More information in the paper :\n
`Zhao, C. and Wang, Z. (2018) GOGO: An improved algorithm to measure the semantic similarity between gene ontology terms. Scientific Reports, 8, 15107; doi:10.1038/s41598-018-33219-y.`
"""

__author__ = "Gérenton Pierre"
__credits__ = ["Gérenton Pierre", "Fabio Zanarello", "Roderic Guigó i Serra"]
__license__ = "CC0 1.0 Universal"
__email__ = "fabio.zanarello@crg.eu"

app = Dash(__name__)

summary = pd.read_csv(args.data + ".summary.tsv", sep = '\t')
data = pd.read_csv(args.data + ".data.tsv",  sep = '\t')


app.layout = html.Div([
    html.H1("Interactive table for intra-gene isoforms diversity"),
    dcc.Tabs(children =[
        dcc.Tab(label="Information", children=[
            html.Br(),
            dcc.Markdown(MARKDOWN_PRESENTATION, mathjax=True)
        ]),
        dcc.Tab(label="Summary", children=[
            html.Br(),
            dash_table.DataTable(data=summary.to_dict('records'),
                                 filter_action='native',
                                 sort_action='native',
                                 columns=[{"name": i, "id": i, "deletable": False, "selectable": False} for i in summary.columns]
                                 ),
        ]),
        dcc.Tab(label="Data", children=[
            html.Br(),
            html.Label('Select columns :'),
            dcc.Dropdown(
                id='column_dropdown',
                options=[{'label': col, 'value': col} for col in data.columns],
                multi=True,
                value=data.columns,
            ),
            html.Br(),
            dash_table.DataTable(id='data_table',
                                 data=data.to_dict('records'),
                                 filter_action='native',
                                 sort_action='native',
                                 columns=[{"name": i, "id": i, "deletable": False, "selectable": False} for i in data.columns],
                                 fixed_rows={'headers': True},
                                 css=[{"selector": ".show-hide", "rule": "display: none"}],
                                 export_format='csv'
                                 ),
        ]),
    ])
])

@callback(
    Output('data_table', 'hidden_columns'),
    Input('column_dropdown', 'value')
)
def update_visible_column_data(selected):
    return list(set(data.columns).difference(set(selected)))

app.run(debug=True)



if __name__=="__main__":
    app.run(debug=True)

