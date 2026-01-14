"""
ProtAtlasViz - Protein Atlas Tissue Expression Visualizer
A Dash application for visualizing tissue-selective gene expression patterns
"""
import dash
from dash import dcc, html, Input, Output, State
import dash_bootstrap_components as dbc
import plotly.graph_objects as go
import plotly.express as px
import pandas as pd
import numpy as np
from data_loader import DataLoader

# Initialize the app
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
app.title = "ProtAtlasViz - Gene Expression Visualizer"

# Initialize data loader
data_loader = DataLoader()
print("Loading expression data...")
data_loader.load_expression_data()
print("Loading histology dictionary...")
data_loader.load_histology_dictionary()
print(f"Data loaded: {len(data_loader.genes)} genes available")

# Layout
app.layout = dbc.Container([
    dbc.Row([
        dbc.Col([
            html.H1("🧬 ProtAtlasViz", className="text-center mb-4 mt-4"),
            html.P("Visualize tissue-selective gene expression patterns from RNA-seq data",
                   className="text-center text-muted mb-4")
        ])
    ]),
    
    dbc.Row([
        dbc.Col([
            dbc.Card([
                dbc.CardBody([
                    html.H5("Gene Selection", className="card-title"),
                    html.P("Search and select up to 10 genes to visualize", 
                           className="text-muted small"),
                    dcc.Dropdown(
                        id='gene-selector',
                        options=[{'label': gene, 'value': gene} for gene in data_loader.genes],
                        multi=True,
                        placeholder="Type to search genes...",
                        searchable=True,
                        maxHeight=300,
                        style={'width': '100%'}
                    ),
                    html.Div(id='gene-count', className="mt-2 small text-muted")
                ])
            ], className="mb-4")
        ], width=12)
    ]),
    
    dbc.Row([
        dbc.Col([
            dbc.Card([
                dbc.CardBody([
                    html.H5("Visualization Options", className="card-title"),
                    dbc.Row([
                        dbc.Col([
                            html.Label("Grouping Level:", className="fw-bold"),
                            dcc.RadioItems(
                                id='grouping-toggle',
                                options=[
                                    {'label': ' Individual Tissues', 'value': 'tissue'},
                                    {'label': ' Organ Groups', 'value': 'organ'}
                                ],
                                value='tissue',
                                inline=True,
                                className="ms-2"
                            )
                        ], md=4),
                        dbc.Col([
                            html.Label("Visualization Type:", className="fw-bold"),
                            dcc.Dropdown(
                                id='viz-type',
                                options=[
                                    {'label': 'Heatmap', 'value': 'heatmap'},
                                    {'label': 'Grouped Bar Chart', 'value': 'bar'},
                                    {'label': 'Box Plot', 'value': 'box'}
                                ],
                                value='heatmap',
                                clearable=False,
                                style={'width': '100%'}
                            )
                        ], md=4),
                        dbc.Col([
                            html.Label("Scale:", className="fw-bold"),
                            dcc.RadioItems(
                                id='log-toggle',
                                options=[
                                    {'label': ' Linear', 'value': 'linear'},
                                    {'label': ' Log', 'value': 'log'}
                                ],
                                value='linear',
                                inline=True,
                                className="ms-2"
                            )
                        ], md=4)
                    ])
                ])
            ], className="mb-4")
        ], width=12)
    ]),
    
    dbc.Row([
        dbc.Col([
            dcc.Loading(
                id="loading",
                type="default",
                children=html.Div(id='visualization-container')
            )
        ], width=12)
    ])
], fluid=True, style={'maxWidth': '1400px'})


def create_heatmap(data: pd.DataFrame, ordered_tissues: list, tissue_to_group: dict, 
                   show_groups: bool = True, log_scale: bool = False) -> go.Figure:
    """Create a heatmap visualization"""
    # Pivot data for heatmap
    pivot_data = data.pivot(index='Gene name', columns='Tissue', values='nTPM')
    
    # Reorder columns to match tissue order
    available_tissues = [t for t in ordered_tissues if t in pivot_data.columns]
    pivot_data = pivot_data[available_tissues]
    
    # Create custom hover text
    hover_text = []
    for gene in pivot_data.index:
        row_text = []
        for tissue in pivot_data.columns:
            value = pivot_data.loc[gene, tissue]
            if pd.isna(value):
                row_text.append(f"Gene: {gene}<br>Tissue: {tissue}<br>nTPM: N/A")
            else:
                group = tissue_to_group.get(tissue, "Other")
                if log_scale:
                    # Show both log and original value
                    original_value = 2**value - 1
                    row_text.append(f"Gene: {gene}<br>Tissue: {tissue}<br>Organ: {group}<br>log2(nTPM+1): {value:.2f}<br>nTPM: {original_value:.1f}")
                else:
                    row_text.append(f"Gene: {gene}<br>Tissue: {tissue}<br>Organ: {group}<br>nTPM: {value:.1f}")
        hover_text.append(row_text)
    
    # Improve color contrast with multi-step colorscale
    colorscale = [
        [0.0, 'rgb(255, 255, 255)'],    # White for zero/low
        [0.1, 'rgb(230, 240, 255)'],    # Very light blue
        [0.3, 'rgb(180, 210, 255)'],    # Light blue
        [0.5, 'rgb(100, 150, 255)'],    # Medium blue
        [0.7, 'rgb(50, 100, 200)'],     # Dark blue
        [1.0, 'rgb(0, 0, 139)']         # Navy blue for maximum
    ]
    
    fig = go.Figure(data=go.Heatmap(
        z=pivot_data.values,
        x=pivot_data.columns,
        y=pivot_data.index,
        colorscale=colorscale,
        hovertemplate='%{customdata}<extra></extra>',
        customdata=hover_text,
        colorbar=dict(title="log2(nTPM+1)" if log_scale else "nTPM")
    ))
    
    # Add group separators if showing individual tissues
    if show_groups:
        current_group = None
        for i, tissue in enumerate(pivot_data.columns):
            group = tissue_to_group.get(tissue, "Other")
            if group != current_group:
                fig.add_vline(x=i-0.5, line_width=2, line_color="white")
                current_group = group
    
    # Calculate dynamic width based on number of tissues (minimum 800, maximum 2000)
    num_tissues = len(pivot_data.columns)
    dynamic_width = max(800, min(2000, num_tissues * 25 + 200))
    
    fig.update_layout(
        title={
            'text': "Gene Expression Heatmap" + (" (Log Scale)" if log_scale else ""),
            'font': {'size': 24, 'color': 'black', 'family': 'Arial'}
        },
        xaxis_title="Tissue",
        yaxis_title="Gene",
        height=max(600, len(pivot_data.index) * 60),
        width=dynamic_width,
        xaxis={
            'tickangle': -45,
            'tickfont': {'size': 14, 'color': 'black', 'family': 'Arial'},
            'titlefont': {'size': 18, 'color': 'black', 'family': 'Arial'}
        },
        yaxis={
            'tickfont': {'size': 14, 'color': 'black', 'family': 'Arial'},
            'titlefont': {'size': 18, 'color': 'black', 'family': 'Arial'}
        },
        font={'family': 'Arial', 'size': 14, 'color': 'black'},
        plot_bgcolor='white',
        paper_bgcolor='white'
    )
    
    return fig


def create_bar_chart(data: pd.DataFrame, ordered_tissues: list, tissue_to_group: dict, 
                     log_scale: bool = False) -> go.Figure:
    """Create a grouped bar chart"""
    # Ensure proper ordering
    data = data.copy()
    data['Tissue'] = pd.Categorical(data['Tissue'], categories=ordered_tissues, ordered=True)
    data = data.sort_values('Tissue')
    
    fig = go.Figure()
    
    for gene in data['Gene name'].unique():
        gene_data = data[data['Gene name'] == gene]
        
        if log_scale:
            hover_text = [
                f"Gene: {gene}<br>Tissue: {tissue}<br>Organ: {tissue_to_group.get(tissue, 'Other')}<br>log2(nTPM+1): {ntpm:.2f}<br>nTPM: {(2**ntpm - 1):.1f}"
                for tissue, ntpm in zip(gene_data['Tissue'], gene_data['nTPM'])
            ]
        else:
            hover_text = [
                f"Gene: {gene}<br>Tissue: {tissue}<br>Organ: {tissue_to_group.get(tissue, 'Other')}<br>nTPM: {ntpm:.1f}"
                for tissue, ntpm in zip(gene_data['Tissue'], gene_data['nTPM'])
            ]
        
        fig.add_trace(go.Bar(
            name=gene,
            x=gene_data['Tissue'],
            y=gene_data['nTPM'],
            hovertemplate='%{customdata}<extra></extra>',
            customdata=hover_text
        ))
    
    yaxis_title = "log2(nTPM+1)" if log_scale else "nTPM"
    
    # Calculate dynamic width based on number of tissues (minimum 800, maximum 2000)
    num_tissues = len(data['Tissue'].unique())
    dynamic_width = max(800, min(2000, num_tissues * 25 + 200))
    
    fig.update_layout(
        title={
            'text': "Gene Expression by Tissue" + (" (Log Scale)" if log_scale else ""),
            'font': {'size': 24, 'color': 'black', 'family': 'Arial'}
        },
        xaxis_title="Tissue",
        yaxis_title=yaxis_title,
        barmode='group',
        height=700,
        width=dynamic_width,
        xaxis={
            'tickangle': -45,
            'tickfont': {'size': 14, 'color': 'black', 'family': 'Arial'},
            'titlefont': {'size': 18, 'color': 'black', 'family': 'Arial'}
        },
        yaxis={
            'tickfont': {'size': 14, 'color': 'black', 'family': 'Arial'},
            'titlefont': {'size': 18, 'color': 'black', 'family': 'Arial'},
            'gridcolor': 'lightgray',
            'gridwidth': 1
        },
        font={'family': 'Arial', 'size': 14, 'color': 'black'},
        plot_bgcolor='white',
        paper_bgcolor='white',
        legend=dict(
            orientation="v",
            yanchor="top",
            y=1,
            xanchor="left",
            x=1.02,
            font={'size': 14, 'color': 'black', 'family': 'Arial'},
            bgcolor='rgba(255,255,255,0.8)',
            bordercolor='black',
            borderwidth=1
        )
    )
    
    return fig


def create_box_plot(data: pd.DataFrame, ordered_tissues: list, tissue_to_group: dict,
                    log_scale: bool = False) -> go.Figure:
    """Create a box plot showing expression distribution"""
    # Ensure proper ordering
    data = data.copy()
    data['Tissue'] = pd.Categorical(data['Tissue'], categories=ordered_tissues, ordered=True)
    data = data.sort_values('Tissue')
    
    fig = go.Figure()
    
    for gene in data['Gene name'].unique():
        gene_data = data[data['Gene name'] == gene]
        
        fig.add_trace(go.Box(
            name=gene,
            x=gene_data['Tissue'],
            y=gene_data['nTPM'],
            boxmean='sd'
        ))
    
    yaxis_title = "log2(nTPM+1)" if log_scale else "nTPM"
    
    # Calculate dynamic width based on number of tissues (minimum 800, maximum 2000)
    num_tissues = len(data['Tissue'].unique())
    dynamic_width = max(800, min(2000, num_tissues * 25 + 200))
    
    fig.update_layout(
        title={
            'text': "Gene Expression Distribution" + (" (Log Scale)" if log_scale else ""),
            'font': {'size': 24, 'color': 'black', 'family': 'Arial'}
        },
        xaxis_title="Tissue",
        yaxis_title=yaxis_title,
        height=700,
        width=dynamic_width,
        xaxis={
            'tickangle': -45,
            'tickfont': {'size': 14, 'color': 'black', 'family': 'Arial'},
            'titlefont': {'size': 18, 'color': 'black', 'family': 'Arial'}
        },
        yaxis={
            'tickfont': {'size': 14, 'color': 'black', 'family': 'Arial'},
            'titlefont': {'size': 18, 'color': 'black', 'family': 'Arial'},
            'gridcolor': 'lightgray',
            'gridwidth': 1
        },
        font={'family': 'Arial', 'size': 14, 'color': 'black'},
        plot_bgcolor='white',
        paper_bgcolor='white',
        legend=dict(
            orientation="v",
            yanchor="top",
            y=1,
            xanchor="left",
            x=1.02,
            font={'size': 14, 'color': 'black', 'family': 'Arial'},
            bgcolor='rgba(255,255,255,0.8)',
            bordercolor='black',
            borderwidth=1
        )
    )
    
    return fig


@app.callback(
    Output('gene-count', 'children'),
    Input('gene-selector', 'value')
)
def update_gene_count(selected_genes):
    """Update the gene count display"""
    if not selected_genes:
        return "No genes selected"
    
    count = len(selected_genes)
    if count > 10:
        return html.Span(f"⚠️ {count}/10 genes selected (maximum 10 allowed)", 
                        style={'color': 'red', 'fontWeight': 'bold'})
    else:
        return f"✓ {count}/10 genes selected"


@app.callback(
    Output('visualization-container', 'children'),
    [Input('gene-selector', 'value'),
     Input('grouping-toggle', 'value'),
     Input('viz-type', 'value'),
     Input('log-toggle', 'value')]
)
def update_visualization(selected_genes, grouping, viz_type, log_scale):
    """Update the visualization based on user selections"""
    if not selected_genes:
        return dbc.Alert(
            "👆 Please select one or more genes to visualize",
            color="info",
            className="text-center"
        )
    
    if len(selected_genes) > 10:
        return dbc.Alert(
            "⚠️ Please select no more than 10 genes",
            color="warning",
            className="text-center"
        )
    
    # Get expression data for selected genes
    expression_data = data_loader.get_expression_for_genes(selected_genes)
    
    if expression_data.empty:
        return dbc.Alert(
            "No expression data found for selected genes",
            color="warning",
            className="text-center"
        )
    
    # Get ordered tissues and grouping information
    ordered_tissues, tissue_to_group = data_loader.get_ordered_tissues_by_group()
    
    # Aggregate by organ group if requested
    if grouping == 'organ':
        expression_data = data_loader.aggregate_by_organ_group(expression_data)
        # Update ordered list for organ groups
        ordered_tissues = sorted(expression_data['Tissue'].unique())
        tissue_to_group = {t: t for t in ordered_tissues}  # Groups map to themselves
    
    # Apply log transformation if requested
    if log_scale == 'log':
        expression_data = expression_data.copy()
        expression_data['nTPM'] = np.log2(expression_data['nTPM'] + 1)  # log2(x+1) to handle zeros
    
    # Create appropriate visualization
    if viz_type == 'heatmap':
        fig = create_heatmap(expression_data, ordered_tissues, tissue_to_group, 
                            show_groups=(grouping == 'tissue'), log_scale=(log_scale == 'log'))
    elif viz_type == 'bar':
        fig = create_bar_chart(expression_data, ordered_tissues, tissue_to_group, 
                              log_scale=(log_scale == 'log'))
    elif viz_type == 'box':
        fig = create_box_plot(expression_data, ordered_tissues, tissue_to_group,
                             log_scale=(log_scale == 'log'))
    else:
        return dbc.Alert("Invalid visualization type", color="danger")
    
    return dcc.Graph(
        figure=fig,
        config={
            'displayModeBar': True,
            'displaylogo': False,
            'toImageButtonOptions': {
                'format': 'png',
                'filename': 'gene_expression_plot',
                'height': 1080,
                'width': None,  # Use figure width
                'scale': 2
            },
            'modeBarButtonsToAdd': ['drawline', 'drawopenpath', 'eraseshape']
        }
    )


if __name__ == '__main__':
    print("\n" + "="*60)
    print("Starting ProtAtlasViz Server")
    print("="*60)
    print("Open your browser and navigate to: http://127.0.0.1:8050/")
    print("="*60 + "\n")
    app.run_server(debug=True, host='127.0.0.1', port=8050)
