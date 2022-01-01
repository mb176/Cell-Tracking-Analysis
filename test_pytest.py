import pytest
from analysis_library import *

def test_read_xml():
    exp = experiment('test')
    exp.read_xml('test_tracks.xml',1)
    assert len(exp.tracks)==2
    assert np.linalg.norm(exp.tracks[1])==np.linalg.norm(np.array([[0,0,0],[1,1,1],[2,2,2]]))
    assert exp.x_max == 3
    assert exp.y_max == 2
    

def test_TAMSD():
    exp = experiment('test')
    exp.read_xml('test_tracks.xml',1)
    #test first track
    assert exp.TAMSD(exp.tracks[0],1) == 1
    #test second track
    assert exp.TAMSD(exp.tracks[1],1) == 2
    #check restriction on delta_t
    with pytest.raises(AssertionError, match='Track duration is shorter than the time interval'):
        exp.TAMSD(exp.tracks[1],2)


def test_average_TAMSD():
    exp = experiment('test')
    exp.read_xml('test_tracks.xml',1)
    #Average of both tracks
    assert exp.average_TAMSD(1,2) == 1.5
    #Only first track 
    assert exp.average_TAMSD(1,4) == 1


#def    test_plot_TAMSD

def test_calculate_tortuosity():
    exp = experiment('test')
    exp.read_xml('test_tracks.xml',1)
    exp.read_xml('test_tracks_1.xml',1) #get extra track 
    exp.calculate_tortuosity()
    assert len(exp.tortuosity) == 3
    assert exp.tortuosity[0] == 0 #linear motion
    assert exp.tortuosity[1] == 0 #linear motion
    assert exp.tortuosity[2] == 1 #returns to starting 

def test_plot_radial_density():
    exp = experiment('test')
    exp.read_xml('test_tracks.xml',1)
    exp.read_xml('test_tracks_1.xml',1)
    histogram, bins = exp.plot_radial_density('axes', 2, 4, no_plot=True)
    assert histogram[0] == 1
    assert histogram[2] == 2

