import json

g_data = {
        "alpha" : 10 ,
        "W_cutofffreq" : 10 ,
        "L_transmissionlen" : 10,
        "prop_speed" : 10
        }

with open('g_data.json', 'w') as f:
    json.dump(g_data, f, indent=2)