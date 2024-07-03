import numpy as np
import copy

def remove_mask_notation(landscape, mask_symbol="-", real_symbol="."):
    """
    Replace mask notation with the normal one, i.e.

    remove_mask_notation(
       ["o#-",
        "..-",
        "*#-"]
       )

    where "-", the mask chracters, are replaced
    by ".", returning:

    ["o#.",
     "...",
     "*#."]
    """
    m = copy.deepcopy(landscape) # the mask

    if isinstance(landscape, str):
        m = m.replace(mask_symbol, real_symbol)
    else:
        for i in range(len(landscape)):
            m[i] = m[i].replace(mask_symbol, real_symbol)
    return m

def get_mask(landscape, mask_symbol="-"):
    """
    returns the mask corresponding to the landscape.
    For instance, 
    
    get_mask(
       ["o#-",
        "..-",
        "*#-"]
       )

    where "-" chracters represent 
    the masked cells, would return:
        
    [0,0,1,
     0,0,1,
     0,0,1]
    
    as a numpy array.
    """
    m = [] # the mask

    if isinstance(landscape, str):
        for i in range(len(landscape)):
            if landscape[i] == mask_symbol:
                m.append(1.)
            else:
                m.append(0.)
    else:
        for i in range(len(landscape)):
            for j in range(len(landscape[i])):
                if landscape[i][j] == mask_symbol:
                    m.append(1.)
                else:
                    m.append(0.)

    return np.array(m)

# tests

assert np.array_equal(get_mask("*o#.--"), np.array([0,0,0,0,1,1]))
assert np.array_equal(get_mask(["*o#",".--"]), np.array([0,0,0,0,1,1]))
assert np.array_equal(get_mask(np.array(["*o#",".--"])), np.array([0,0,0,0,1,1]))

assert remove_mask_notation("*o#.--") == "*o#..."
assert remove_mask_notation(["*o#",".--"]) == ["*o#","..."]
assert np.array_equal(remove_mask_notation(np.array(["*o#",".--"])), np.array(["*o#","..."]))

assert np.array_equal(get_mask("*o#.$$", "$"), np.array([0,0,0,0,1,1]))
assert np.array_equal(get_mask(["*o#",".$$"], "$"), np.array([0,0,0,0,1,1]))
assert np.array_equal(get_mask(np.array(["*o#",".$$"]), "$"), np.array([0,0,0,0,1,1]))

assert remove_mask_notation("*o#.$$", "$", "*") == "*o#.**"
assert remove_mask_notation(["*o#",".$$"], "$", "*") == ["*o#",".**"]
assert np.array_equal(remove_mask_notation(np.array(["*o#",".$$"]), "$", "*"), np.array(["*o#",".**"]))

