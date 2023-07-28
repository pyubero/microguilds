#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 12:15:36 2023

@author: juanrivassantisteban
"""

import pandas as pd
import plotly.offline as pyo
import plotly.graph_objs as go
import numpy as np
import Levenshtein

# Read data from the TSV file
try:
    df = pd.read_csv("kvalues_n_uptake_Species_GTDB.tsv", sep="\t")
except FileNotFoundError:
    print("Error: File not found. Make sure the file 'kvalues_n_uptake_Species_GTDB.tsv' is in the correct path.")
    exit()

# Convert "Cluster" and "Context" columns to string data type
df['Cluster'] = df['Cluster'].astype(str)
df['Context'] = df['Context'].astype(str)

# Convert the "k-value" column to numeric, setting any invalid values to NaN
df['k-value'] = pd.to_numeric(df['k-value'], errors='coerce')

# Drop rows with NaN values (caused by invalid data in numeric conversion)
df = df.dropna()

# Take the logarithm of the "k-value" column (base 10 in this example)
df['log_k-value'] = df['k-value'].apply(lambda x: np.log10(x))

# Scale the 'log_k-value' to adjust marker size (you can adjust the scale factor as needed)
scale_factor = 6
df['marker_size'] = df['log_k-value'].apply(lambda x: np.abs(x) * scale_factor)

# Create a dictionary to map Taxa to colors
taxa_colors = {
    'Pseudomonas': 'red',
    'Marinobacter': 'yellow',
    'Alcanivorax': 'green',
    'Punicei': 'pink',
    'Methylobac': 'black',
    'Cyclocl': 'black',
    'Salipig': 'blue',
    'Bacteroid': 'red',
    'Gamma': 'yellow',
    'Nitrosp': 'green',
    'Cyano': 'green',
    'Beta': 'pink',
    'Methylobac': 'black',
    'Cyclocl': 'black',
    'Alpha': 'blue'
    
}

# Set default color to blue for all other Taxa not in the dictionary
df['marker_color'] = df['Taxon'].apply(lambda x: next((color for taxon, color in taxa_colors.items() if taxon in x), 'grey'))

# Function to calculate the Levenshtein similarity score between two strings
def calculate_similarity(name, target):
    max_distance = max(len(name), len(target))
    levenshtein_distance = Levenshtein.distance(name, target)
    similarity_score = 1 - levenshtein_distance / max_distance
    return similarity_score

# Target name for similarity comparison
target_taxon = 'Pseudomonas'  # You can change this to any other taxon for comparison

# Filter DataFrame to include only rows where "Context" contains "MP"
df = df[df['Context'].str.contains('MP')]

# Calculate similarity scores for each taxa name in comparison to the target_taxon
df['similarity'] = df['Taxon'].apply(lambda x: calculate_similarity(x, target_taxon))

# Sort the DataFrame by similarity scores in descending order (higher similarity first)
df = df.sort_values(by='similarity', ascending=False)

# Create 3D scatter plot using Plotly
fig = go.Figure()

# Add 3D scatter plot with hover information
fig.add_trace(go.Scatter3d(x=df['Context'], y=df['Cluster'], z=df['Taxon'], mode='markers',
                           marker=dict(size=df['marker_size'], color=df['marker_color'],
                                       colorscale='Viridis', colorbar=dict(title='Taxa')),
                           text=df['Taxon']))  # Set 'text' property for hover information

# Get unique context values
unique_contexts = df['Context'].unique()

# Create planes for each unique context
fig.update_layout(scene=dict(zaxis=dict(showticklabels=False)))



x_coords = np.arange(100, len(unique_contexts))  # Generate equally spaced x-coordinates

for x_coord, context_value in zip(x_coords, unique_contexts):
    df_context = df[df['Context'] == context_value]
    plane = go.Mesh3d(x=[x_coord]*4, y=df_context['Cluster'].unique(), z=df_context['Taxon'].unique(),
                      i=[0, 0, 0, 0], j=[1, 1, 2, 2], k=[2, 3, 3, 2],
                      opacity=1, color='rgba(700, 200, 200, 1)')
    fig.add_trace(plane)

fig.update_layout(scene=dict(xaxis=dict(title=''), yaxis=dict(title='Cluster'), zaxis=dict(title='Taxon')),
                  margin=dict(l=0, r=0, b=0, t=0))



# Display the plot in the Spyder IDE
pyo.plot(fig)
