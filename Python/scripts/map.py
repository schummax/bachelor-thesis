import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

x = 14

# Increase the image size
fig, ax = plt.subplots(figsize=(8, 6), edgecolor="w", dpi=300)

# Create the main map (m)
m = Basemap(
    projection="merc",
    resolution="f",
    llcrnrlat=59,
    urcrnrlat=73,
    llcrnrlon=120,
    urcrnrlon=170,
    ax=ax,
)
m.drawcountries()
# m.fillcontinents()  # color="lightgray"
# m.etopo(scale=0.5, alpha=0.5)
m.shadedrelief(scale=0.5, alpha=0.5)

plt.xlabel(labelpad=20, xlabel="Longitude")
plt.ylabel(labelpad=20, ylabel="Latitude")

# Draw latitude and longitude lines for the main map
parallels = range(40, 81, 10)
meridians = range(20, 181, 20)
m.drawparallels(parallels, labels=[1, 0, 0, 0], fontsize=10)
m.drawmeridians(meridians, labels=[0, 0, 0, 1], fontsize=10)

# Provided coordinates
locations = [
    (68.37016, 161.41555),
    (61.75967, 130.47438),
    (61.76490, 130.46503),
    (61.76086, 130.47466),
]

# Plot the provided coordinates on the main map
for i, (lat, lon) in enumerate(locations):
    if i == 0:
        x, y = m(lon, lat)
        m.plot(x, y, "bo", markersize=3)
        plt.text(x, y, f"Loc {i}", fontsize=x, ha="left", va="bottom")
        continue
    x, y = m(lon, lat)
    m.plot(x, y, "bo", markersize=3)


# Create an inset plot for the zoomed-in map (g)
axins = inset_axes(ax, width="50%", height="50%", loc="upper left", borderpad=3)
g = Basemap(
    projection="merc",
    resolution="f",
    llcrnrlat=61.75,
    urcrnrlat=61.77,
    llcrnrlon=130.45,
    urcrnrlon=130.49,
    ax=axins,
)
g.drawcountries()
g.fillcontinents(color="#cad7d4")  # color="lightgray"
# g.shadedrelief(scale=0.5, alpha=0.5)


# Draw latitude and longitude lines for the inset map
inset_parallels = np.arange(61.75, 61.77, 0.01)
inset_meridians = np.arange(130.47, 130.48, 0.01)
g.drawparallels(inset_parallels, labels=[1, 0, 0, 0], fontsize=7)
g.drawmeridians(inset_meridians, labels=[0, 0, 0, 1], fontsize=7)

# Plot the provided coordinates on the zoomed-in map
for i, (lat, lon) in enumerate(locations[1:]):
    if i == 0:
        x, y = g(lon, lat)
        g.plot(x, y, "bo", markersize=3)
        axins.text(x, y, f"Loc {i+1}", fontsize=x, ha="center", va="top")
        continue
    x, y = g(lon, lat)
    g.plot(x, y, "bo", markersize=3)
    axins.text(x, y, f"Loc {i+1}", fontsize=x, ha="center", va="bottom")

# Draw lines from location 3 to the bottom corners of the inset map
loc3_x, loc3_y = m(locations[3][1], locations[3][0])
llcrnrx, llcrnry = m(127.5, 65.63)
urcrnrx, urcrnry = m(144.3, 65.63)

m.plot([loc3_x, llcrnrx], [loc3_y, llcrnry], color="black", linewidth=0.5)
m.plot([loc3_x, urcrnrx], [loc3_y, urcrnry], color="black", linewidth=0.5)

plt.savefig("../output/map.png", bbox_inches="tight")
