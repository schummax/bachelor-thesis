import altair as alt
import pandas as pd

# Read in the data
data = pd.read_csv('AncientMetagenomeDir_Rigou.tsv', sep='\t')


map = alt.topo_feature('https://github.com/zarkzork/russia-topojson/blob/master/russia.json', 'states')

alt.Chart(map).mark_geoshape(
    fill='lightgray',
    stroke='white'
)

# # Create a map of the sample locations
# points = alt.Chart(data).mark_circle(size=60).encode(
#     y='longitude:Q',
#     x='latitude:Q',
#     tooltip=['sample_name', 'latitude', 'longitude']
# )

# chart = background + points
# chart.save('sample_locations_map.html')