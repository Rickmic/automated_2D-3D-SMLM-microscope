"""
Some helpful functions for the processing and analysis of 2D and 3D localization data. These functions are
compatible with the LocData class of locan [1]. The AIM drift correction function was original publisched in [2].
The implementation is based on the implementation in Picasso https://github.com/jungmannlab/picasso.
(MIT license, Copyright (c) 2016 Jungmann Lab, MPI of Biochemistry) [3].

References
----------
..  [1] S. Doose, LOCAN: a python library for
    analyzing single-molecule localization
    microscopy data, Bioinformatics, 2022

    [2] Hongqiang Ma, Maomao Chen, Phuong Nguyen, Yang Liu,
    Toward drift-free high-throughput nanoscopy through
    adaptive intersection maximization, Science Advances, 2024

    [3] J. Schnitzbauer*, M.T. Strauss*, T. Schlichthaerle,
    F. Schueder, R. Jungmann, Super-Resolution Microscopy
    with DNA-PAINT, Nature Protocols, 2017
"""


from __future__ import annotations
import cv2
import math
import numpy
import numpy as np
import trackpy as tp
import pandas as pd
import locan as lc
import logging
import copy
import time
import open3d as o3d
from sklearn.neighbors import NearestNeighbors
import matplotlib.pyplot as plt
from collections.abc import Callable, Sequence
from typing import TYPE_CHECKING, Any, Literal
from locan.data import LocData
from locan import Bins, histogram
from locan.configuration import COLORMAP_DEFAULTS
from locan.visualize.colormap import get_colormap
from locan.visualize.transform import adjust_contrast
from concurrent.futures import ThreadPoolExecutor
from scipy.interpolate import InterpolatedUnivariateSpline
from shapely.geometry import LineString
from shapely import BufferCapStyle, BufferJoinStyle
from tqdm import tqdm


if TYPE_CHECKING:
    import boost_histogram as bh
    import matplotlib as mpl

    from locan.visualize.transform import Trafo

logger = logging.getLogger(__name__)




def track_beads(
        locdata: LocData,
        search_range: int,
        memory: int,
        bead_locs: int,
        three_dimension: bool = False,
) -> LocData:
    """
    Track beads and sample with trackpy trajectories:

    Parameters
    ----------
    locdata
        locan locdata object containing beads and sample localizations
    search_range
        the maximum distance features can move between frames (trackpy.link_df)
    memory
        the maximum number of frames during which a feature can vanish,
        then reappear nearby, and be considered the same particle. (trackpy.link_df)
    bead_locs
        number of localizations a bead should have at least
    three_dimension
        if true localisations are tracked in 3D, default False

    Returns
    -------
    locdata_beads
        A new LocData instance with localization data of detected beads
    locdata_sampleFloat
        A new LocData instance with localization data of the sample without beads
    """


    if three_dimension:
        trajectories = tp.link_df(locdata.data, search_range=search_range, memory=memory, pos_columns=['position_x', 'position_y', 'position_z'], t_column='frame')
    else:
        trajectories = tp.link_df(locdata.data, search_range=search_range, memory=memory,
                                  pos_columns=['position_x', 'position_y'], t_column='frame')
    particle_abundance = trajectories['particle'].value_counts().to_frame() # list with all particle counts
    bead_df = particle_abundance.query(f"count > {bead_locs}") # DataFrame with all beads
    bead_list = list(bead_df.index.values)
    print(f"{len(bead_list)} beads found with more than {bead_locs} localizations")
    mask = trajectories['particle'].isin(bead_list) # generate boolean mask for beads
    beads = trajectories[mask]
    sample = trajectories[~mask]
    locdata_beads = lc.LocData.from_dataframe(dataframe=beads)
    locdata_sample = lc.LocData.from_dataframe(dataframe=sample)

    return locdata_beads, locdata_sample


def extract_beads(
        bead_file: LocData,
        three_dimension: bool = False,
) -> pd.DataFrame:
    """
    Extract mean x and y positions of tracked beads.

    Parameters
    ----------
    bead_file
        Localization data of tracked beads with particle column containing the particle number
    three_dimensions
        if true localisations are tracked in 3D, default False

    Returns
    -------
    pandas DataFrame with particle number, mean coordinates
    and the size of the trajectories in x and y
    """
    particles = pd.DataFrame(bead_file.data["particle"].unique(), columns=["particle_number"])

    x_mean = []
    x_size = []
    y_mean = []
    y_size = []
    z_mean = []
    z_size = []

    if three_dimension:
        for particle in particles["particle_number"]:
            bead = bead_file.data.loc[bead_file.data['particle'] == particle]
            x_mean.append(bead["position_x"].mean())
            x_size.append(bead["position_x"].max() - bead["position_x"].min()) # calculate trajectorie size
            y_mean.append(bead["position_y"].mean())
            y_size.append(bead["position_y"].max() - bead["position_y"].min()) # calculate trajectorie size
            z_mean.append(bead["position_z"].mean())
            z_size.append(bead["position_z"].max() - bead["position_z"].min())  # calculate trajectorie size
        particles["x_mean"] = x_mean
        particles["x_size"] = x_size
        particles["y_mean"] = y_mean
        particles["y_size"] = y_size
        particles["z_mean"] = z_mean
        particles["z_size"] = z_size
    else:
        for particle in particles["particle_number"]:
            bead = bead_file.data.loc[bead_file.data['particle'] == particle]
            x_mean.append(bead["position_x"].mean())
            x_size.append(bead["position_x"].max() - bead["position_x"].min()) # calculate trajectorie size
            y_mean.append(bead["position_y"].mean())
            y_size.append(bead["position_y"].max() - bead["position_y"].min()) # calculate trajectorie size
        particles["x_mean"] = x_mean
        particles["x_size"] = x_size
        particles["y_mean"] = y_mean
        particles["y_size"] = y_size

    return particles

def discard_beads(
        bead_position: pd.DataFrame,
        min_distance: int,
        three_dimension: bool = False,
) -> pd.DataFrame:
    """
    Discard beads, that are closer than min distance to each other.

    Parameters
    ----------
    bead_position
        pandas DataFrame with particle number, mean coordinates
        and the size of the trajectories in x and y
    min_distance
        Minimum distance in nanometer two beads should be seperated at least
    three_dimensions
        if true localisations are tracked in 3D, default False
    Returns
    -------
    pandas DataFrame with additional column discard with boolean values
    """

    bead_position["discard"] = [False] * len(bead_position)
    if three_dimension:
        for index in bead_position.index:
            # calculate coordinates for actual bead:
            index_coordinates = [bead_position["x_mean"][index], bead_position["y_mean"][index], bead_position["z_mean"][index]]
            #calculate coordinates for all beads, that should be compared:
            for compare in range((index + 1), (len(bead_position))):
                compare_coordinates = [bead_position["x_mean"][compare], bead_position["y_mean"][compare], bead_position["z_mean"][compare]]
                distance = math.dist(index_coordinates, compare_coordinates)
                if distance <= min_distance:
                    bead_position.at[compare, "discard"] = True
                    bead_position.at[index, "discard"] = True
                    print(f"Beads {index} and {compare} are only {distance} nm seperated and will be discarded")
    else:
        for index in bead_position.index:
            # calculate coordinates for actual bead:
            index_coordinates = [bead_position["x_mean"][index], bead_position["y_mean"][index]]
            #calculate coordinates for all beads, that should be compared:
            for compare in range((index + 1), (len(bead_position))):
                compare_coordinates = [bead_position["x_mean"][compare], bead_position["y_mean"][compare]]
                distance = math.dist(index_coordinates, compare_coordinates)
                if distance <= min_distance:
                    bead_position.at[compare, "discard"] = True
                    bead_position.at[index, "discard"] = True
                    print(f"Beads {index} and {compare} are only {distance} nm seperated and will be discarded")

    return bead_position



def affine_transf_matrix(
        bead_positions_1: pd.DataFrame,
        bead_positions_2: pd.DataFrame,
        max_distance: int,
        three_dimension: bool = False,
) -> np.ndarray:
    """
    Calculate matrix for affine transformation based on fiducial coordinates. bead_positions_1 and bead_positions_2 should
    contain only beads ,that are seperated by more than max_distance.

    Parameters
    ----------
    bead_positions_1
        pandas DataFrame with at least columns x_mean and y_mean representing the X and Y coordinates of the fiducials from color Channel 1.
    bead_positions_2
        pandas DataFrame with at least columns x_mean and y_mean representing the X and Y coordinates of the fiducials from color Channel 2.
    max_distance
        Maximal distance between two coordinates from position 1 and 2, at which the coordinates should be grouped together.
    three_dimensions
        if true localisations are tracked in 3D, default False

    Returns
    -------
    numpy ndarray with 2D affine transformation matrix (2x3)
    """
    if three_dimension:
        ## calculate transform coordinates:
        used_beads = 0
        discard_beads = 0
        transform_coordinates_1 = []
        transform_coordinates_2 = []
        for index in bead_positions_1.index:
            # calculate coordinates from bead_position_1:
            pos1_coordinates = [bead_positions_1["x_mean"][index], bead_positions_1["y_mean"][index], bead_positions_1["z_mean"][index]]
            for compare in bead_positions_2.index:
                pos2_coordinates = [bead_positions_2["x_mean"][compare], bead_positions_2["y_mean"][compare], bead_positions_2["z_mean"][compare]]
                color_distance = math.dist(pos1_coordinates, pos2_coordinates)
                if color_distance <= max_distance: # beads are corresponding
                    if (not bead_positions_1["discard"][index]) and (not bead_positions_2["discard"][compare]):
                        transform_coordinates_1.append(pos1_coordinates)
                        transform_coordinates_2.append(pos2_coordinates)
                        used_beads += 1
                else:
                    discard_beads += 1
        print(
            f"{used_beads} beads were used to create the affine transformation matrix and {discard_beads} beads were discarded")
        ## generate affine transformation matrix:
        matrix = cv2.estimateAffine3D(np.array(transform_coordinates_1), np.array(transform_coordinates_2))[1]
    else:
        ## calculate transform coordinates:
        used_beads = 0
        discard_beads = 0
        transform_coordinates_1 = []
        transform_coordinates_2 = []
        for index in bead_positions_1.index:
            # calculate coordinates from bead_position_1:
            pos1_coordinates = [bead_positions_1["x_mean"][index], bead_positions_1["y_mean"][index]]
            for compare in bead_positions_2.index:
                pos2_coordinates = [bead_positions_2["x_mean"][compare], bead_positions_2["y_mean"][compare]]
                color_distance = math.dist(pos1_coordinates, pos2_coordinates)
                if color_distance <= max_distance:  # beads are corresponding
                    if (not bead_positions_1["discard"][index]) and (not bead_positions_2["discard"][compare]):
                        transform_coordinates_1.append(pos1_coordinates)
                        transform_coordinates_2.append(pos2_coordinates)
                        used_beads += 1
                else:
                    discard_beads += 1
        print(f"{used_beads} beads were used to create the affine transformation matrix and {discard_beads} beads were discarded")
        ## generate affine transformation matrix:
        matrix = cv2.estimateAffine2D(np.array(transform_coordinates_1), np.array(transform_coordinates_2))[0]
    return matrix, used_beads








def render_2d_np(
    locdata: LocData,
    loc_properties: list[str] | None = None,
    other_property: str | None = None,
    bins: Bins | bh.axis.Axis | bh.axis.AxesTuple | None = None,
    n_bins: int | Sequence[int] | None = None,
    bin_size: float | Sequence[float] | Sequence[Sequence[float]] | None = 10,
    bin_edges: Sequence[float] | Sequence[Sequence[float]] | None = None,
    bin_range: (
        tuple[float, float]
        | Sequence[float]
        | Sequence[Sequence[float]]
        | Literal["zero", "link"]
        | None
    ) = None,
    rescale: int | str | Trafo | Callable[..., Any] | bool | None = None,
) -> numpy.ndarray:
    """
    Render localization data into a 2D numpy array by binning x,y-coordinates into
    regular bins. Based on: render_2d_mpl from locan

    Parameters
    ----------
    locdata
        Localization data.
    loc_properties
        Localization properties to be grouped into bins.
        If None The coordinate_values of `locdata` are used.
    other_property
        Localization property (columns in `locdata.data`) that is averaged
        in each pixel.
        If None, localization counts are shown.
    bins
        The bin specification as defined in :class:`Bins`
    bin_edges
        Bin edges for all or each dimension
        with shape (dimension, n_bin_edges).
    bin_range
        Minimum and maximum edge for all or each dimensions
        with shape (2,) or (dimension, 2).
        If None (min, max) ranges are determined from data and returned;
        if 'zero' (0, max) ranges with max determined from data are returned.
        if 'link' (min_all, max_all) ranges with min and max determined from
        all combined data are returned.
    n_bins
        The number of bins for all or each dimension.
        5 yields 5 bins in all dimensions.
        (2, 5) yields 2 bins for one dimension and 5 for the other dimension.
    bin_size
        The size of bins for all or each bin and for all or each dimension
        with shape (dimension,) or (dimension, n_bins).
        5 would describe bin_size of 5 for all bins in all dimensions.
        ((2, 5),) yield bins of size (2, 5) for one dimension.
        (2, 5) yields bins of size 2 for one dimension and 5 for the other
        dimension.
        ((2, 5), (1, 3)) yields bins of size (2, 5) for one dimension and
        (1, 3) for the other dimension.
        To specify arbitrary sequence of `bin_size` use `bin_edges` instead.
    rescale
        Transformation as defined in :class:`locan.Trafo` or by
        transformation function.
        For None or False no rescaling occurs.
        Legacy behavior:
        For tuple with upper and lower bounds provided in percent,
        rescale intensity values to be within percentile of max and min
        intensities.
        For 'equal' intensity values are rescaled by histogram equalization.

    Returns
    -------
    numpy.ndarray
        ndarray with the image.
    """

    # create warning if single localization
    if len(locdata) == 1:
        logger.warning("Locdata carries a single localization.")

    data, bins, labels = histogram(
        locdata=locdata,
        loc_properties=loc_properties,
        other_property=other_property,
        bins=bins,
        n_bins=n_bins,
        bin_size=bin_size,
        bin_edges=bin_edges,
        bin_range=bin_range,
    )
    data = adjust_contrast(data, rescale)


    return data

def save_2d_png(
    locdata: LocData,
    path: str,
    cmap: str | None = None,
    loc_properties: list[str] | None = None,
    other_property: str | None = None,
    bins: Bins | bh.axis.Axis | bh.axis.AxesTuple | None = None,
    n_bins: int | Sequence[int] | None = None,
    bin_size: float | Sequence[float] | Sequence[Sequence[float]] | None = 10,
    bin_edges: Sequence[float] | Sequence[Sequence[float]] | None = None,
    bin_range: (
        tuple[float, float]
        | Sequence[float]
        | Sequence[Sequence[float]]
        | Literal["zero", "link"]
        | None
    ) = None,
    rescale: int | str | Trafo | Callable[..., Any] | bool | None = None,
):
    """
    Render localization data into a 2D numpy array by binning x,y-coordinates into
    regular bins. Based on: render_2d_mpl from locan

    Parameters
    ----------
    locdata
        Localization data.
    path
        filepath for png file
    loc_properties
        Localization properties to be grouped into bins.
        If None The coordinate_values of `locdata` are used.
    other_property
        Localization property (columns in `locdata.data`) that is averaged
        in each pixel.
        If None, localization counts are shown.
    bins
        The bin specification as defined in :class:`Bins`
    bin_edges
        Bin edges for all or each dimension
        with shape (dimension, n_bin_edges).
    bin_range
        Minimum and maximum edge for all or each dimensions
        with shape (2,) or (dimension, 2).
        If None (min, max) ranges are determined from data and returned;
        if 'zero' (0, max) ranges with max determined from data are returned.
        if 'link' (min_all, max_all) ranges with min and max determined from
        all combined data are returned.
    n_bins
        The number of bins for all or each dimension.
        5 yields 5 bins in all dimensions.
        (2, 5) yields 2 bins for one dimension and 5 for the other dimension.
    bin_size
        The size of bins for all or each bin and for all or each dimension
        with shape (dimension,) or (dimension, n_bins).
        5 would describe bin_size of 5 for all bins in all dimensions.
        ((2, 5),) yield bins of size (2, 5) for one dimension.
        (2, 5) yields bins of size 2 for one dimension and 5 for the other
        dimension.
        ((2, 5), (1, 3)) yields bins of size (2, 5) for one dimension and
        (1, 3) for the other dimension.
        To specify arbitrary sequence of `bin_size` use `bin_edges` instead.
    rescale
        Transformation as defined in :class:`locan.Trafo` or by
        transformation function.
        For None or False no rescaling occurs.
        Legacy behavior:
        For tuple with upper and lower bounds provided in percent,
        rescale intensity values to be within percentile of max and min
        intensities.
        For 'equal' intensity values are rescaled by histogram equalization.

    Returns
    -------
    """

    # create warning if single localization
    if len(locdata) == 1:
        logger.warning("Locdata carries a single localization.")

    data, bins, labels = histogram(
        locdata=locdata,
        loc_properties=loc_properties,
        other_property=other_property,
        bins=bins,
        n_bins=n_bins,
        bin_size=bin_size,
        bin_edges=bin_edges,
        bin_range=bin_range,
    )
    data = adjust_contrast(data, rescale)

    if not cmap:
        lc_cmap = COLORMAP_DEFAULTS["CONTINUOUS"]
        cmap = get_colormap(colormap=lc_cmap).matplotlib

    plt.imsave(path, data, cmap=cmap)


def intersect1d(a, b):
    """Slightly faster implementation of _np.intersect1d without
    unnecessary checks, etc.

    Finds the indices of common elements in two 1D arrays (a and b).
    Both a and b are assumed to be sorted and contain only unique
    values.

    Parameters
    ----------
    a : _np.array
        1D array of integers.
    b : _np.array
        1D array of integers.

    Returns
    -------
    a_indices : _np.array
        Indices of common elements in a.
    b_indices : _np.array
        Indices of common elements in b.
    """

    aux = np.concatenate((a, b))
    aux_sort_indices = np.argsort(aux, kind='mergesort')
    aux = aux[aux_sort_indices]

    mask = aux[1:] == aux[:-1]
    a_indices = aux_sort_indices[:-1][mask]
    b_indices = aux_sort_indices[1:][mask] - a.size

    return a_indices, b_indices


def count_intersections(l0_coords, l0_counts, l1_coords, l1_counts):
    """Counts the number of intersected localizations between the two
    datasets. We assume that the intersection distance is 1 and since
    the coordinates are expressed in the units of intersection distance,
    we require the coordinates to be exactly the same to count as
    intersection. Also, coordinates are converted to 1D arrays
    (x + y * width).

    Parameters
    ----------
    l0_coords : _np.array
        Unique coordinates of the reference localizations.
    l0_counts : _np.array
        Counts of the unique values of reference localizations.
    l1_coords : _np.array
        Unique coordinates of the target localizations.
    l1_counts : _np.array
        Counts of the unique values of target localizations.

    Returns
    -------
    n_intersections : int
        Number of intersections.
    """

    # indices of common elements
    idx0, idx1 = intersect1d(l0_coords, l1_coords)
    # extract the counts of these elements
    l0_counts_subset = l0_counts[idx0]
    l1_counts_subset = l1_counts[idx1]
    # for each overlapping coordinate, take the minimum count from l0
    # and l1, sum up across all overlapping coordinates
    n_intersections = np.sum(
        np.minimum(l0_counts_subset, l1_counts_subset)
    )
    return n_intersections


def run_intersections(
        l0_coords, l0_counts, l1_coords, l1_counts, shifts_xy, box
):
    """Run intersection counting across the local search region. Returns
    the 2D array with number of intersections across the local search
    region.

    Parameters
    ----------
    l0_coords : _np.array
        Unique coordinates of the reference localizations.
    l0_counts : _np.array
        Counts of the reference localizations.
    l1_coords : _np.array
        Unique coordinates of the target localizations.
    l1_counts : _np.array
        Counts of the target localizations.
    shifts_xy : _np.array
        1D array with x and y shifts.
    box : int
        Side length of the local search region.

    Returns
    -------
    roi_cc : _np.2darray
        2D array with number of intersections across the local search
        region.
    """

    # create the 2D array with shifts
    roi_cc = np.zeros(shifts_xy.shape, dtype=np.int32)
    # shift target coordinates
    l1_coords_shifted = l1_coords[:, np.newaxis] + shifts_xy
    # go through each element in the local search region
    for i in range(len(shifts_xy)):
        n_intersections = count_intersections(
            l0_coords, l0_counts, l1_coords_shifted[:, i], l1_counts
        )
        roi_cc[i] = n_intersections
    return roi_cc.reshape(box, box)


def run_intersections_multithread(
        l0_coords, l0_counts, l1_coords, l1_counts, shifts_xy, box
):
    """Run intersection counting across the local search region. Returns
    the 2D array with number of intersections across the local search
    region. Uses multithreading.

    Parameters
    ----------
    l0_coords : _np.array
        Unique coordinates of the reference localizations.
    l0_counts : _np.array
        Counts of the reference localizations.
    l1_coords : _np.array
        Unique coordinates of the target localizations.
    l1_counts : _np.array
        Counts of the target localizations.
    shifts_xy : _np.array
        1D array with x and y shifts.
    box : int
        Side length of the local search region.

    Returns
    -------
    roi_cc : _np.2darray
        2D array with number of intersections across the local search
        region.
    """

    # shift target coordinates
    l1_coords_shifted = l1_coords[:, np.newaxis] + shifts_xy
    # run multiple threads
    n_workers = len(shifts_xy)
    executor = ThreadPoolExecutor(n_workers)
    f = [
        executor.submit(
            count_intersections,
            l0_coords,
            l0_counts,
            l1_coords_shifted[:, i],
            l1_counts,
        )
        for i in range(len(shifts_xy))
    ]
    executor.shutdown(wait=True)
    if box == 1:  # z intersection only, for z undrifting
        roi_cc = np.array([_.result() for _ in f])
    else:  # 2D intersection
        roi_cc = np.array([_.result() for _ in f]).reshape(box, box)
    return roi_cc


def point_intersect_2d(
        l0_coords, l0_counts, x1, y1, intersect_d, width_units, shifts_xy, box
):
    """Converts target coordinates into a 1D array in units of
    intersect_d and counts the number of intersections in the local
    search region.

    Parameters
    ----------
    l0_coords : _np.array
        Unique values of the reference localizations.
    l0_counts : _np.array
        Counts of the unique values of reference localizations.
    x1 : _np.array
        x-coordinates of the target (currently undrifted) localizations.
    y1 : _np.array
        y-coordinates of the target (currently undrifted) localizations.
    intersect_d : float
        Intersect distance in camera pixels.
    width_units : int
        Width of the camera image in units of intersect_d.
    shifts_xy : _np.array
        1D array with x and y shifts.
    box : int
        Final side length of the local search region.

    Returns
    -------
    roi_cc : _np.2darray
        2D array with numbers of intersections in the local search
        region.
    """

    # convert target coordinates to a 1D array in intersect_d units
    x1_units = np.round(x1 / intersect_d)
    y1_units = np.round(y1 / intersect_d)
    l1 = np.int32(x1_units + y1_units * width_units)  # 1d list
    # get unique values and counts of the target localizations
    l1_coords, l1_counts = np.unique(l1, return_counts=True)
    # run the intersections counting
    roi_cc = run_intersections_multithread(
        l0_coords, l0_counts, l1_coords, l1_counts, shifts_xy, box
    )
    return roi_cc


def point_intersect_3d(
        l0_coords, l0_counts, x1, y1, z1,
        intersect_d, width_units, height_units, shifts_z,
):
    """Converts target coordinates into a 1D array in units of
    intersect_d and counts the number of intersections in the local
    search region.

    Parameters
    ----------
    l0_coords : _np.array
        Unique values of the reference localizations.
    l0_counts : _np.array
        Counts of the unique values of reference localizations.
    x1 : _np.array
        x-coordinates of the target (currently undrifted) localizations.
    y1 : _np.array
        y-coordinates of the target (currently undrifted) localizations.
    z1 : _np.array
        z-coordinates of the target (currently undrifted) localizations.
    intersect_d : float
        Intersect distance in camera pixels.
    width_units : int
        Width of the camera image in units of intersect_d.
    height_units : int
        Height of the camera image in units of intersect_d.
    shifts_z : _np.array
        1D array with z shifts.

    Returns
    -------
    roi_cc : _np.2darray
        2D array with numbers of intersections in the local search
        region.
    """

    # convert target coordinates to a 1D array in intersect_d units
    x1_units = np.round(x1 / intersect_d)
    y1_units = np.round(y1 / intersect_d)
    z1_units = np.round(z1 / intersect_d)
    l1 = np.int32(
        x1_units
        + y1_units * width_units
        + z1_units * width_units * height_units
    )  # 1d list
    # get unique values and counts of the target localizations
    l1_coords, l1_counts = np.unique(l1, return_counts=True)
    # run the intersections counting
    roi_cc = run_intersections_multithread(
        l0_coords, l0_counts, l1_coords, l1_counts, shifts_z, 1
    )
    return roi_cc


def get_fft_peak(roi_cc, roi_size):
    """Estimate the precise sub-pixel position of the peak of roi_cc
    with FFT.

    Parameters
    ----------
    roi_cc : _np.2darray
        2D array with numbers of intersections in the local search region.
    roi_size : int
        Size of the local search region.

    Returns
    -------
    px : float
        Estimated x-coordinate of the peak.
    py : float
        Estimated y-coordinate of the peak.
    """

    fft_values = np.fft.fft2(roi_cc.T)
    ang_x = np.angle(fft_values[0, 1])
    ang_x = ang_x - 2 * np.pi * (ang_x > 0)  # normalize
    px = (
            np.abs(ang_x) / (2 * np.pi / roi_cc.shape[0])
            - (roi_cc.shape[0] - 1) / 2
    )  # peak in x
    px *= roi_size / roi_cc.shape[0]  # convert to intersect_d units
    ang_y = np.angle(fft_values[1, 0])
    ang_y = ang_y - 2 * np.pi * (ang_y > 0)  # normalize
    py = (
            np.abs(ang_y) / (2 * np.pi / roi_cc.shape[1])
            - (roi_cc.shape[1] - 1) / 2
    )  # peak in y
    py *= roi_size / roi_cc.shape[1]  # convert to intersect_d units
    return px, py



def get_fft_peak_z(roi_cc, roi_size):
    """Estimate the precise sub-pixel position of the peak of 1D roi_cc.

    Parameters
    ----------
    roi_cc : _np.array
        1D array with numbers of intersections in the local search
        region.
    roi_size : int
        Size of the local search region.

    Returns
    -------
    pz : float
        Estimated z-coordinate of the peak.
    """

    fft_values = np.fft.fft(roi_cc)
    ang_z = np.angle(fft_values[1])
    ang_z = ang_z - 2 * np.pi * (ang_z > 0)  # normalize
    pz = (
            np.abs(ang_z) / (2 * np.pi / roi_cc.size)
            - (roi_cc.size - 1) / 2
    )  # peak in z
    pz *= roi_size / roi_cc.size  # convert to intersect_d units
    return pz



def intersection_max(
        x, y, ref_x, ref_y,
        frame, seg_bounds, intersect_d, roi_r, width,
        aim_round=1, progress=None,
):
    """Maximize intersection (undrift) for 2D localizations.

    Parameters
    ----------
    x : _np.array
        x-coordinates of the localizations.
    y : _np.array
        y-coordinates of the localizations.
    ref_x_list : _np.array
        x-coordinates of the reference localizations.
    ref_y_list : _np.array
        y-coordinates of the reference localizations.
    frame : _np.array
        Frame indices of localizations.
    seg_bounds : np.array
        Frame indices of the segmentation bounds. Defines temporal
        intervals used to estimate drift.
    intersect_d : float
        Intersect distance in camera pixels.
    roi_r : float
        Radius of the local search region in camera pixels. Should be
        higher than the maximum expected drift within one segment.
    width : int
        Width of the camera image in camera pixels.
    aim_round : {1, 2}
        Round of AIM algorithm. The first round uses the first interval
        as reference, the second round uses the entire dataset as
        reference. The impact is that in the second round, the first
        interval is also undrifted.
    progress : picasso.lib.ProgressDialog (default=None)
        Progress dialog. If None, progress is displayed with tqdm.

    Returns
    -------
    x_pdc : _np.array
        Undrifted x-coordinates.
    y_pdc : _np.array
        Undrifted y-coordinates.
    drift_x : _np.array
        Drift in x-direction.
    drift_y : _np.array
        Drift in y-direction.
    """

    assert aim_round in [1, 2], "aim_round must be 1 or 2."

    # number of segments
    n_segments = len(seg_bounds) - 1
    rel_drift_x = 0  # adaptive drift (updated at each interval)
    rel_drift_y = 0

    # drift in x and y
    drift_x = np.zeros(n_segments)
    drift_y = np.zeros(n_segments)

    # find shifts for the local search region (in units of intersect_d)
    roi_units = int(np.ceil(roi_r / intersect_d))
    steps = np.arange(-roi_units, roi_units + 1, 1)
    box = len(steps)
    shifts_xy = np.zeros((box, box), dtype=np.int32)
    width_units = width / intersect_d
    for i, shift_x in enumerate(steps):
        for j, shift_y in enumerate(steps):
            shifts_xy[i, j] = shift_x + shift_y * width_units
    shifts_xy = shifts_xy.reshape(box ** 2)

    # convert reference to a 1D array in units of intersect_d and find
    # unique values and counts
    x0_units = np.round(ref_x / intersect_d)
    y0_units = np.round(ref_y / intersect_d)
    l0 =np.int32(x0_units + y0_units * width_units)  # 1d list
    l0_coords, l0_counts = np.unique(l0, return_counts=True)

    # initialize progress such that if GUI is used, tqdm is omitted
    start_idx = 1 if aim_round == 1 else 0
    if progress is not None:
        iterator = range(start_idx, n_segments)
    else:
        iterator = tqdm(
            range(start_idx, n_segments),
            desc=f"Undrifting ({aim_round}/2)",
            unit="segment",
        )

        # run across each segment
    for s in iterator:
        # get the target localizations within the current segment
        min_frame_idx = frame > seg_bounds[s]
        max_frame_idx = frame <= seg_bounds[s + 1]
        x1 = x[min_frame_idx & max_frame_idx]
        y1 = y[min_frame_idx & max_frame_idx]

        # skip if no reference localizations
        if len(x1) == 0:
            drift_x[s] = drift_x[s - 1]
            drift_y[s] = drift_y[s - 1]
            continue

        # undrifting from the previous round
        x1 += rel_drift_x
        y1 += rel_drift_y

        # count the number of intersected localizations
        roi_cc = point_intersect_2d(
            l0_coords, l0_counts, x1, y1,
            intersect_d, width_units, shifts_xy, box,
        )

        # estimate the precise sub-pixel position of the peak of roi_cc
        # with FFT
        px, py = get_fft_peak(roi_cc, 2 * roi_r)

        # update the relative drift reference for the subsequent
        # segmented subset (interval) and save the drifts
        rel_drift_x += px
        rel_drift_y += py
        drift_x[s] = -rel_drift_x
        drift_y[s] = -rel_drift_y

        # update progress
        if progress is not None:
            progress.set_value(s)
        else:
            iterator.update(s - iterator.n)

    # interpolate the drifts (cubic spline) for all frames
    t = (seg_bounds[1:] + seg_bounds[:-1]) / 2
    drift_x_pol = InterpolatedUnivariateSpline(t, drift_x, k=3)
    drift_y_pol = InterpolatedUnivariateSpline(t, drift_y, k=3)
    t_inter = np.arange(seg_bounds[-1]) + 1
    drift_x = drift_x_pol(t_inter)
    drift_y = drift_y_pol(t_inter)

    # undrift the localizations
    x_pdc = x - drift_x[frame - 1]
    y_pdc = y - drift_y[frame - 1]

    return x_pdc, y_pdc, drift_x, drift_y


def intersection_max_z(
        x, y, z, ref_x, ref_y, ref_z,
        frame, seg_bounds, intersect_d, roi_r, width, height, pixelsize,
        aim_round=1, progress=None,
):
    """Maximize intersection (undrift) for 3D localizations. Assumes
    that x and y coordinates were already undrifted. x and y are in
    units of camera pixels, z is in nm.

    See intersection_max for more details."""

    # convert z to camera pixels
    z = z.copy() / pixelsize
    ref_z = ref_z.copy() / pixelsize  # TODO: remember to convert back to nm, also for drift_z

    # number of segments
    n_segments = len(seg_bounds) - 1
    rel_drift_z = 0  # adaptive drift (updated at each interval)

    # drift in z
    drift_z = np.zeros(n_segments)

    # find shifts for the local search region (in units of intersect_d)
    roi_units = int(np.ceil(roi_r / intersect_d))
    steps = np.arange(-roi_units, roi_units + 1, 1)
    box = len(steps)
    width_units = width / intersect_d
    height_units = height / intersect_d
    shifts_z = steps.astype(np.int32) * width_units * height_units

    # convert reference to a 1D array in units of intersect_d and find
    # unique values and counts
    x0_units = np.round(ref_x / intersect_d)
    y0_units = np.round(ref_y / intersect_d)
    z0_units = np.round(ref_z / intersect_d)
    l0 = np.int32(
        x0_units
        + y0_units * width_units
        + z0_units * width_units * height_units
    )  # 1d list
    l0_coords, l0_counts = np.unique(l0, return_counts=True)

    # initialize progress such that if GUI is used, tqdm is omitted
    start_idx = 1 if aim_round == 1 else 0
    if progress is not None:
        iterator = range(start_idx, n_segments)
    else:
        iterator = tqdm(
            range(start_idx, n_segments),
            desc=f"Undrifting z ({aim_round}/2)",
            unit="segment",
        )

    # run across each segment
    for s in iterator:
        # get the target localizations within the current segment
        min_frame_idx = frame > seg_bounds[s]
        max_frame_idx = frame <= seg_bounds[s + 1]
        x1 = x[min_frame_idx & max_frame_idx]
        y1 = y[min_frame_idx & max_frame_idx]
        z1 = z[min_frame_idx & max_frame_idx]

        # skip if no reference localizations
        if len(x1) == 0:
            drift_z[s] = drift_z[s - 1]
            continue

        # undrifting from the previous round
        z1 += rel_drift_z

        # count the number of intersected localizations
        roi_cc = point_intersect_3d(
            l0_coords, l0_counts, x1, y1, z1,
            intersect_d, width_units, height_units, shifts_z,
        )

        # estimate the precise sub-pixel position of the peak of roi_cc
        # with FFT
        pz = get_fft_peak_z(roi_cc, 2 * roi_r)

        # update the relative drift reference for the subsequent
        # segmented subset (interval) and save the drifts
        rel_drift_z += pz
        drift_z[s] = -rel_drift_z

        # update progress
        if progress is not None:
            progress.set_value(s)
        else:
            iterator.update(s - iterator.n)

    # interpolate the drifts (cubic spline) for all frames
    t = (seg_bounds[1:] + seg_bounds[:-1]) / 2
    drift_z_pol = InterpolatedUnivariateSpline(t, drift_z, k=3)
    t_inter = np.arange(seg_bounds[-1]) + 1
    drift_z = drift_z_pol(t_inter)

    # undrift the localizations
    z_pdc = z - drift_z[frame - 1]

    # convert back to nm
    z_pdc *= pixelsize
    drift_z *= pixelsize

    return z_pdc, drift_z



def estimate_drift_aim(
        locdata: LocData,
        segmentation: int = 100,
        pixelsize: float = 111,
        width: int = 256,
        height: int = 256,
        n_frames: int = 5_000,
        first_frame: int = 0,
        intersect_d: float = 20 / 111,
        roi_r: float = 60 / 111,
        intersect_d_z: float = 40 / 111,
        roi_r_z: float = 60 / 111,
        progress=None,
        thirdD: bool = True,
) -> (LocData, np.rec.array):
    """Apply AIM undrifting to the localizations.

    Parameters
    ----------
    locs
        Localizations list to be undrifted.
    intersect_d
        Intersect distance in camera pixels.
    segmentation
        Time interval for drift tracking, unit: frames.
    roi_r
        Radius of the local search region in camera pixels. Should be
        larger than the  maximum expected drift within segmentation.
    progress : picasso.lib.ProgressDialog (default=None)
        Progress dialog. If None, progress is displayed with into the
        console.

    Returns
    -------
    locs
        Undrifted localizations.
    drift
        Drift in x and y directions (and z if applicable).
    """

    if np.isnan(width * height ):
        raise KeyError(
            "Insufficient metadata available. Please specify 'Width', 'Height',"
            " 'Frames' and 'Pixelsize' in the metadata .yaml."
        )
    if np.isnan(pixelsize * n_frames):
        raise KeyError(
            "Insufficient metadata available. Please specify 'Width', 'Height',"
            " 'Frames' and 'Pixelsize' in the metadata .yaml."
        )
    # x an y coordinates in camera pixels:
    locdata.data["position_x"] = locdata.data["position_x"]/pixelsize
    locdata.data["position_y"] = locdata.data["position_y"]/pixelsize

    # frames should start at 1
    frame = locdata.data["frame"] + 1 - locdata.data["frame"].min()

    # find the segmentation bounds (temporal intervals)
    seg_bounds = np.concatenate((
        np.arange(0, n_frames, segmentation), [n_frames]
    ))

    # get the reference localizations (first interval)
    ref_x = locdata.data["position_x"][frame <= segmentation]
    ref_y = locdata.data["position_y"][frame <= segmentation]

    ### RUN AIM TWICE ###
    # the first run is with the first interval as reference
    x_pdc, y_pdc, drift_x1, drift_y1 = intersection_max(
        locdata.data.position_x, locdata.data.position_y, ref_x, ref_y,
        frame, seg_bounds, intersect_d, roi_r, width,
        aim_round=1, progress=progress,
    )
    # the second run is with the entire dataset as reference
    if progress is not None:
        progress.zero_progress(description="Undrifting by AIM (2/2)")
    x_pdc, y_pdc, drift_x2, drift_y2 = intersection_max(
        x_pdc, y_pdc, x_pdc, y_pdc,
        frame, seg_bounds, intersect_d, roi_r, width,
        aim_round=2, progress=progress,
    )

    # add the drifts together from the two rounds
    drift_x = drift_x1 + drift_x2
    drift_y = drift_y1 + drift_y2

    # shift the drifts by the mean value
    #shift_x = np.mean(drift_x)
    #shift_y = np.mean(drift_y)
    #drift_x -= shift_x
    #drift_y -= shift_y
    #x_pdc += shift_x
    #y_pdc += shift_y

    # combine to Picasso format
    pos_frame = np.arange(first_frame, first_frame+ n_frames)
    drift_x_px = drift_x* pixelsize
    drift_y_px = drift_y* pixelsize
    drift = np.rec.array((drift_x_px, drift_y_px, pos_frame), dtype=[("x", "f"), ("y", "f"), ("frame", "f")])

    # 3D undrifting
    if hasattr(locdata.data, "position_z") and thirdD:
        if progress is not None:
            progress.zero_progress(description="Undrifting z (1/2)")
        ref_x = x_pdc[frame <= segmentation]
        ref_y = y_pdc[frame <= segmentation]
        ref_z = locdata.data.position_z[frame <= segmentation]
        z_pdc, drift_z1 = intersection_max_z(
            x_pdc, y_pdc, locdata.data.position_z, ref_x, ref_y, ref_z,
            frame, seg_bounds, intersect_d_z, roi_r_z, width, height, pixelsize,
            aim_round=1, progress=progress,
        )
        if progress is not None:
            progress.zero_progress(description="Undrifting z (2/2)")
        z_pdc, drift_z2 = intersection_max_z(
            x_pdc, y_pdc, z_pdc, x_pdc, y_pdc, z_pdc,
            frame, seg_bounds, intersect_d_z, roi_r_z, width, height, pixelsize,
            aim_round=2, progress=progress,
        )
        drift_z = drift_z1 + drift_z2
        shift_z = np.mean(drift_z)
        drift_z -= shift_z
        z_pdc += shift_z
        drift = np.rec.array(
            (drift_x, drift_y, drift_z),
            dtype=[("x", "f"), ("y", "f"), ("z", "f")]
        )

    # apply the drift to localizations
    locdata.data["position_x"] = x_pdc * pixelsize
    locdata.data["position_y"] = y_pdc * pixelsize
    if hasattr(locdata.data, "position_z") and thirdD:
        locdata.data["position_z"] = z_pdc



    if progress is not None:
        progress.close()

    return (locdata, drift)


def get_coord_columns(df):
    """Find coordinate columns (3D or 2D) in a LocData dataframe.

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe from a LocData object containing localization data.

    Returns
    -------
    tuple of str
        Tuple with names of the x, y, z coordinate columns. If only 2D
        coordinates are present, the third entry (z) is returned as
        ``None``.

    Raises
    ------
    ValueError
        If no suitable coordinate columns can be found.
    """
    candidates_3d = [
        ("x", "y", "z"),
        ("position_x", "position_y", "position_z"),
        ("X", "Y", "Z"),
    ]
    for cols in candidates_3d:
        if all(c in df.columns for c in cols):
            return cols

    candidates_2d = [
        ("x", "y"),
        ("position_x", "position_y"),
        ("X", "Y"),
    ]
    for cols in candidates_2d:
        if all(c in df.columns for c in cols):
            return cols + (None,)

    raise ValueError("No suitable coordinate columns found in dataframe.")


def locdata_to_array(ld):
    """Convert a LocData object to a NumPy array of coordinates.

    Parameters
    ----------
    ld : locan.LocData
        LocData object containing localizations.

    Returns
    -------
    arr : numpy.ndarray
        Array of shape (N, 3) containing x, y, z coordinates. If the
        input is 2D, the z column is filled with zeros.
    coord_cols : tuple of str
        Tuple with the names of the x, y, z coordinate columns used.
        If only 2D information was present, the z entry is ``None``.
    """
    df = ld.data
    xcol, ycol, zcol = get_coord_columns(df)

    if zcol is None:
        arr = np.stack(
            [
                df[xcol].to_numpy(),
                df[ycol].to_numpy(),
                np.zeros(len(df), dtype=np.float32),
            ],
            axis=1,
        ).astype(np.float32)
    else:
        arr = df[[xcol, ycol, zcol]].to_numpy(dtype=np.float32)

    return arr, (xcol, ycol, zcol)


def rmse_nn(A, B):
    """Compute RMSE to nearest neighbors between two point sets.

    Parameters
    ----------
    A : numpy.ndarray
        Array of shape (N, D) representing source points.
    B : numpy.ndarray
        Array of shape (M, D) representing target points.

    Returns
    -------
    float
        Root mean square error of distances from each point in ``A``
        to its nearest neighbor in ``B``.
    """
    nn = NearestNeighbors(n_neighbors=1).fit(B)
    dists, _ = nn.kneighbors(A)
    return float(np.sqrt(np.mean(dists**2)))


def align_locdata_ransac_icp(
    moving_ld: lc.LocData,
    fixed_ld: lc.LocData,
    localization_sigma: float = 20.0,
    icp_dist_factors=(5.0, 3.0, 1.5),
    ransac_max_iter: int = 100000,
    ransac_confidence: float = 0.999,
    align_z: bool = True,
    verbose: bool = True,
):
    """Align localizations via RANSAC and multi-scale ICP.

    The moving LocData is globally aligned to the fixed LocData using
    FPFH-based RANSAC for coarse registration followed by multi-scale
    ICP refinement on the full point clouds.

    Parameters
    ----------
    moving_ld : locan.LocData
        LocData object containing the localizations to be aligned
        (moving / source).
    fixed_ld : locan.LocData
        LocData object containing the reference localizations
        (fixed / target).
    localization_sigma : float, optional
        Approximate localization precision in the same units as the
        coordinates (e.g. nm). Used to set multi-scale ICP distance
        thresholds. Default is 20.0.
    icp_dist_factors : tuple of float, optional
        Multiplicative factors applied to ``localization_sigma`` to
        define the maximum correspondence distance at each ICP stage
        (coarse → fine). Default is (5.0, 3.0, 1.5).
    ransac_max_iter : int, optional
        Maximum number of RANSAC iterations for global registration.
        Default is 100000.
    ransac_confidence : float, optional
        Confidence level for RANSAC convergence. Default is 0.999.
    align_z : bool, optional
        If True, perform full 3D alignment (x, y, z). If False,
        the registration is performed in 2D (x, y) by projecting
        all points to z=0 for RANSAC/ICP, and only x and y
        coordinates of ``moving_ld`` are updated. The original
        z coordinates are preserved. Default is True.
    verbose : bool, optional
        If True, print progress and diagnostics to the console.
        Default is True.

    Returns
    -------
    aligned_ld : locan.LocData
        New LocData object containing the aligned localizations of
        ``moving_ld`` in the coordinate system of ``fixed_ld``.
    T_icp : numpy.ndarray
        Final 4×4 homogeneous transformation matrix estimated from
        the registration (applied to points with z=0 if
        ``align_z=False``).
    """
    # --------------------------------------------------
    # 1) LocData → NumPy (full point clouds)
    # --------------------------------------------------

    pc1, cols1 = locdata_to_array(moving_ld)
    pc2, cols2 = locdata_to_array(fixed_ld)
    x1, y1, z1 = cols1

    if verbose:
        print("pc1 (moving) shape:", pc1.shape)
        print("pc2 (fixed) shape:", pc2.shape)

    # For registration, optionally project to 2D (z=0)
    if align_z:
        reg_pc1 = pc1.copy()
        reg_pc2 = pc2.copy()
    else:
        reg_pc1 = pc1.copy()
        reg_pc2 = pc2.copy()
        reg_pc1[:, 2] = 0.0
        reg_pc2[:, 2] = 0.0
        if verbose:
            print("align_z=False: using z=0 projection for registration.")

    pcd_src = o3d.geometry.PointCloud()
    pcd_tgt = o3d.geometry.PointCloud()
    pcd_src.points = o3d.utility.Vector3dVector(reg_pc1)
    pcd_tgt.points = o3d.utility.Vector3dVector(reg_pc2)

    # Bounding box of target for parameter scaling (use registration coords)
    bbox_min = reg_pc2.min(axis=0)
    bbox_max = reg_pc2.max(axis=0)
    bbox_diag = np.linalg.norm(bbox_max - bbox_min)
    if bbox_diag <= 0:
        bbox_diag = 1.0

    if verbose:
        print(f"Bounding box diagonal (target, registration coords): {bbox_diag:.3f}")

    # Just for sanity
    c1 = reg_pc1.mean(axis=0)
    c2 = reg_pc2.mean(axis=0)
    if verbose:
        print("Centroid distance before registration (registration coords):",
              np.linalg.norm(c1 - c2))

    # --------------------------------------------------
    # 2) Global registration with RANSAC (Open3D)
    # --------------------------------------------------

    if verbose:
        print("\n=== Global registration: RANSAC (Open3D) ===")

    # 2a) Estimate normals + FPFH features on full point clouds
    radius_normal = bbox_diag / 50.0
    radius_feature = bbox_diag / 10.0

    pcd_src.estimate_normals(
        o3d.geometry.KDTreeSearchParamHybrid(radius=radius_normal, max_nn=30)
    )
    pcd_tgt.estimate_normals(
        o3d.geometry.KDTreeSearchParamHybrid(radius=radius_normal, max_nn=30)
    )

    src_fpfh = o3d.pipelines.registration.compute_fpfh_feature(
        pcd_src,
        o3d.geometry.KDTreeSearchParamHybrid(radius=radius_feature, max_nn=100),
    )
    tgt_fpfh = o3d.pipelines.registration.compute_fpfh_feature(
        pcd_tgt,
        o3d.geometry.KDTreeSearchParamHybrid(radius=radius_feature, max_nn=100),
    )

    # Slightly looser threshold for RANSAC (global)
    # You can also use: ransac_distance_thresh = 10 * localization_sigma
    ransac_distance_thresh = bbox_diag / 50.0

    t0 = time.perf_counter()
    ransac_result = o3d.pipelines.registration.registration_ransac_based_on_feature_matching(
        pcd_src,
        pcd_tgt,
        src_fpfh,
        tgt_fpfh,
        mutual_filter=True,
        max_correspondence_distance=ransac_distance_thresh,
        estimation_method=o3d.pipelines.registration.TransformationEstimationPointToPoint(),
        ransac_n=4,
        checkers=[
            o3d.pipelines.registration.CorrespondenceCheckerBasedOnEdgeLength(0.9),
            o3d.pipelines.registration.CorrespondenceCheckerBasedOnDistance(ransac_distance_thresh),
        ],
        criteria=o3d.pipelines.registration.RANSACConvergenceCriteria(
            max_iteration=ransac_max_iter,
            confidence=ransac_confidence,
        ),
    )
    t1 = time.perf_counter()

    if verbose:
        print(f"RANSAC runtime: {t1 - t0:.3f} s")
        print("RANSAC fitness:", ransac_result.fitness,
              "inlier RMSE:", ransac_result.inlier_rmse)
        print("RANSAC transformation:\n", ransac_result.transformation)

    # Sanity check: reject obviously bad RANSAC transforms
    T_ransac = ransac_result.transformation
    ransac_translation = np.linalg.norm(T_ransac[:3, 3])

    if (ransac_result.fitness < 0.02) or (ransac_translation > bbox_diag) or (
        ransac_result.inlier_rmse > 5 * localization_sigma
    ):
        if verbose:
            print("RANSAC looks unreliable (low fitness, large translation, or high inlier RMSE).")
            print("Using identity as initial transform for ICP.")
        init_transform = np.eye(4)
    else:
        init_transform = T_ransac

    # --------------------------------------------------
    # 3) Multi-scale ICP refinement
    # --------------------------------------------------

    if verbose:
        print("\n=== Multi-scale ICP refinement ===")
    icp_max_dists = [f * localization_sigma for f in icp_dist_factors]
    if verbose:
        print("ICP max correspondence distances (coarse → fine):", icp_max_dists)

    current_transform = init_transform

    for stage, max_dist in enumerate(icp_max_dists, start=1):
        if verbose:
            print(f"\n--- ICP stage {stage} with max_dist = {max_dist:.2f} ---")
        t0 = time.perf_counter()
        icp_result = o3d.pipelines.registration.registration_icp(
            pcd_src,
            pcd_tgt,
            max_dist,
            current_transform,
            o3d.pipelines.registration.TransformationEstimationPointToPoint(),
            o3d.pipelines.registration.ICPConvergenceCriteria(
                max_iteration=80
            )
        )
        t1 = time.perf_counter()
        current_transform = icp_result.transformation

        if verbose:
            print(f"ICP stage {stage} runtime: {t1 - t0:.3f} s")
            print(f"ICP stage {stage} fitness: {icp_result.fitness:.4f}, "
                  f"inlier RMSE: {icp_result.inlier_rmse:.4f}")

    T_icp = current_transform
    if verbose:
        print("\nFinal ICP transformation:\n", T_icp)

    # --------------------------------------------------
    # 4) Apply final transform to original coordinates
    # --------------------------------------------------

    if align_z:
        # Full 3D transform
        pc1_h = np.hstack([pc1, np.ones((pc1.shape[0], 1), dtype=pc1.dtype)])
        pc1_icp = (T_icp @ pc1_h.T).T[:, :3]
        # 3D RMSE
        rmse_icp = rmse_nn(pc1_icp, pc2)
    else:
        # Only align x,y using transform estimated on z=0 projection
        pc1_for_tf = pc1.copy()
        pc1_for_tf[:, 2] = 0.0
        pc1_h = np.hstack([pc1_for_tf, np.ones((pc1_for_tf.shape[0], 1), dtype=pc1_for_tf.dtype)])
        pc1_tf_plane = (T_icp @ pc1_h.T).T[:, :3]

        # Start from original coords, replace x,y by transformed x,y
        pc1_icp = pc1.copy()
        pc1_icp[:, 0:2] = pc1_tf_plane[:, 0:2]  # new x,y
        # z remains unchanged

        # 2D RMSE (xy only)
        rmse_icp = rmse_nn(pc1_icp[:, :2], pc2[:, :2])

    if verbose:
        print(f"\nFinal RANSAC + multi-scale ICP RMSE(src->tgt): {rmse_icp:.3f}")
        if align_z:
            c1_icp = pc1_icp.mean(axis=0)
            print("Centroid distance after ICP (3D):",
                  np.linalg.norm(c1_icp - c2))
        else:
            c1_icp_xy = pc1_icp[:, :2].mean(axis=0)
            c2_xy = pc2[:, :2].mean(axis=0)
            print("Centroid distance after ICP (2D, xy):",
                  np.linalg.norm(c1_icp_xy - c2_xy))

    # --------------------------------------------------
    # 5) New LocData object with transformed coordinates
    # --------------------------------------------------

    df1 = moving_ld.data.copy()
    df1_icp = df1.copy()
    df1_icp[x1] = pc1_icp[:, 0]
    df1_icp[y1] = pc1_icp[:, 1]
    if z1 is not None and align_z:
        df1_icp[z1] = pc1_icp[:, 2]

    aligned_ld = lc.LocData(
        dataframe=df1_icp,
        meta=copy.deepcopy(moving_ld.meta)
    )

    if verbose:
        print("\n=== Done ===")
        print("New aligned LocData object created.")
        print(f"RANSAC + multi-scale ICP RMSE: {rmse_icp:.3f}")
        if not align_z:
            print("Note: alignment was performed in 2D (x, y only); z was kept unchanged.")

    return aligned_ld, T_icp



def centerline_to_polygon(coords, thickness, cap_style="flat", join_style="round"):
    """
    Convert a napari polyline to a polygon of given thickness.

    Parameters
    ----------
    coords : (N, 2) or (N, 3) array
        Polyline coordinates from a napari Shapes layer.
        If 3D, buffering is done in the last two dimensions and the first
        dimension (e.g. z) is kept constant.
    thickness : float
        Total polygon thickness in the same units as the coordinates.
    cap_style : {"round", "flat", "square"}
        End style of the thickened line.
    join_style : {"round", "mitre", "bevel"}
        Corner style.

    Returns
    -------
    polygon : (M, 2) or (M, 3) array
        Polygon coordinates usable in napari.
    """
    coords = np.asarray(coords, dtype=float)

    if coords.ndim != 2 or coords.shape[0] < 2:
        raise ValueError("coords must have shape (N, 2) or (N, 3) with N >= 2")

    if coords.shape[1] < 2:
        raise ValueError("coords must contain at least 2 spatial dimensions")

    # Use the last two dimensions for 2D buffering (typically y, x)
    line_2d = coords[:, -2:]

    cap_map = {
        "round": BufferCapStyle.round,
        "flat": BufferCapStyle.flat,
        "square": BufferCapStyle.square,
    }
    join_map = {
        "round": BufferJoinStyle.round,
        "mitre": BufferJoinStyle.mitre,
        "bevel": BufferJoinStyle.bevel,
    }

    poly = LineString(line_2d).buffer(
        thickness / 2,
        cap_style=cap_map[cap_style],
        join_style=join_map[join_style],
    )

    # If buffering creates multiple polygons, keep the largest one
    if poly.geom_type == "MultiPolygon":
        poly = max(poly.geoms, key=lambda g: g.area)

    polygon_2d = np.asarray(poly.exterior.coords)

    # If original coordinates are 3D, keep the first dimension fixed
    if coords.shape[1] > 2:
        fixed_dims = coords[0, :-2]
        fixed_block = np.repeat(fixed_dims[None, :], len(polygon_2d), axis=0)
        polygon = np.concatenate([fixed_block, polygon_2d], axis=1)
    else:
        polygon = polygon_2d

    return polygon


def line_to_area(line, z_range):
    """
    Calculate contour in 3D from a given line with a certain range in z:

    Parameters
    ----------
    line
        numpy.ndarray with two columns representing x and y coordinates of the line
    z_range
        int representing the z range of the line in nm

    Returns
    -------
    area
        Float representing the area of the line in nm^2
    """
    first_coord = line[0,:]
    sec_coord = line[1:,:]
    distance = 0
    for line in sec_coord:
        distance += math.dist(first_coord, line)
        first_coord = line
    area = distance*z_range
    return area



def clustert_3d(label,
               condition,
               roi_files,
               line_files,
               eps,
               min_samples,
               nanodomain_threshold = 30
               ):
    """
    Perform a DBSCAN in 3D of a given list of rois and calculate the number of clusters, cluster density
    and percentage of large clusters:

    Parameters
    ----------
    label
        str representing the name of the whole experiment
    condition
        str representing the condition of the files in the plots later on
    roi_files
        list containing the file path to all rois (YAML files).
    line_files
        list containing the file paths to all ROI lines (.npy) (corresponding to each ROI).
    eps
        The maximum distance between two samples for them to be considered as
        in the same neighborhood. (Parameter passed to locan cluster_dbscan)
    min_samples
        The number of samples in a neighborhood for a point to be considered
        as a core point. This includes the point itself. (Parameter passed to locan cluster_dbscan)
    nanodomain_threshold
        Threshold in clusters at which the cluster is considered as nanodomain.

    Returns
    -------
    cluster_number
        pd.DataFrame containing the number of clusters for each ROI.
    cluster_area
        pd.DataFrame containing the cluster density for each ROI.
    large_cluster
        pd.DataFrame containing the percentage of large clusters for each ROI.
    pdf_dataframe
        pd.DataFrame containing the cluster size for all clusters.
    """
    cluster_number_list = []
    cluster_area_list = []
    large_cluster_list = []
    pdf_list = []

    selector = list(range(0, len(roi_files)))
    print('Selector:')
    print(selector)
    print('---------------------------------------------')

    for position in selector:
        # load SMLM data
        print(f"Processing file: {roi_files[position]}")
        print(f"with line: {line_files[position]}")
        print('---------------------------------------------')
        roi = lc.Roi.from_yaml(roi_files[position])
        roi_locs = roi.locdata()

        # load center line
        c_line = np.load(line_files[position])

        # cluster SMLM data
        roi_noise, roi_clust = lc.cluster_dbscan(roi_locs, eps=eps, min_samples=min_samples)
        cluster_number_pos = roi_clust.meta.element_count
        cluster_number_list.append([label, condition, cluster_number_pos])
        #pdf_list.extend(roi_clust.data['localization_count'].to_list())
        # Store cluster size together with the selector/ROI index
        for cluster_size in roi_clust.data['localization_count']:
            pdf_list.append([position, cluster_size])

        # calculate area
        roi_area = line_to_area(c_line, 4000)
        cluster_area_pos = cluster_number_pos / (roi_area / 1000_000)
        cluster_area_list.append([label, condition, cluster_area_pos])


        # calculate large cluster:
        large_clust = lc.select_by_condition(roi_clust, condition=f'localization_count > {nanodomain_threshold}')
        large_clust.reduce()
        small_clust = lc.select_by_condition(roi_clust, condition=f'localization_count < {nanodomain_threshold}')
        small_clust.reduce()
        large_clust = (large_clust.data.shape[0] /(small_clust.data.shape[0] + large_clust.data.shape[0]))*100
        large_cluster_list.append([label, condition, large_clust])

    cluster_number = pd.DataFrame(cluster_number_list, columns=['label','condition', 'cluster'],)
    cluster_area = pd.DataFrame(cluster_area_list, columns=['label','condition', 'cluster / µm²'],)
    large_cluster = pd.DataFrame(large_cluster_list, columns=['label','condition', 'nanodomains [%]'],)
    #pdf_dataframe = pd.DataFrame(pdf_list, columns=['cluster size'],)
    pdf_dataframe = pd.DataFrame(pdf_list, columns=['selector', 'cluster size'])

    return cluster_number, cluster_area, large_cluster, pdf_dataframe
