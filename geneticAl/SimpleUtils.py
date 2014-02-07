try: # Works for Python 2.6+
    import json as simplejson
except ImportError: # Otherwise you need to install simplejson
    import simplejson

def besttrans(etens, etoscs, return_osc=False):
    """Return the scaling factor and energy of the best transition"""
    idx = 0
    found = -1
    max = 0
    for idx in range(0, len(etens)):
        if etoscs[idx] >= 1.0 and found==-1:
            found = idx
        if etoscs[idx] > etoscs[max]:
            max = idx

    scale = 1.0
    chosenidx = found
    if found == -1:
        # Use the strongest oscillator, and scale by the
        # strength
        chosenidx = max
        scale = etoscs[chosenidx]

    if not return_osc:
        return scale, etens[chosenidx]
    else:
        return scale, etens[chosenidx], etoscs[chosenidx]

def getHplusBG(json):
    """Get the HOMO and lowest energy significant transition
    given the JSON from the sqlite database"""
    homo, lumo, etens, etoscs = simplejson.loads(json)
    scale, trans = besttrans(etens, etoscs)
    return homo, trans

