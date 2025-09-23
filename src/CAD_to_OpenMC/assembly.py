import cadquery as cq
import numpy as np
import OCP
import pathlib as pl
from typing import List, Optional, Union

from itertools import zip_longest

import re
import os
import math
import sys
from cadquery import BoundBox
from unidecode import unidecode
from datetime import datetime

from pymoab import core, types

import CAD_to_OpenMC.assemblymesher as am
from CAD_to_OpenMC.datadirectory import mesher_datadir
from CAD_to_OpenMC.check_step import has_degenerate_toroids

try:
    import gmsh
    nogmsh = False
except ImportError as e:
    print(f"Warning: import gmsh failed ({e})- material tag-list, must be supplied.")
    nogmsh = True

mesher_config = {
    # general opts
    "default": False,
    "vetoed": None,
    # stl-backend opts
    "tolerance": 0.1,
    "angular_tolerance": 0.2,
    # gmsh-backend opts
    "min_mesh_size": 0.1,
    "max_mesh_size": 10,
    "curve_samples": 20,
    "mesh_algorithm": 1,
    "threads": 1,
    "radial_threshold": 0,
    "refine": 0,
    "verbose": 0,
}


# these are dummies that we still need to have defined
def _replace(filename, string1, string2):
    pass


class Entity:
    """
    Container class for Parts / Volumes

    This class is a container to allow iterating through the geometry model, storing links to
    data as we go along.
    Each Entity is a glue-object that connects a step part with a h5m-volume, including positioning,
    material tags, and links to surface mesh.


    Parameters
    ----------
    solid : cadquery:Solid
        The cadquery Solid object corresponding to the step-part
    idx : int
        The index number of the part/solid.
    tag : string
        Material tag of the part/solid
    """

    def __init__(self, solid=None, idx=0, tag="vacuum"):
        if solid is not None:
            self.solid = solid
            self.idx = idx
            self.tag = tag
            # extract some parameters from the solid for later.
            self.bb = solid.BoundingBox()
            self.center = solid.Center()
            self.volume = solid.Volume()
            # When availalble, the bounding box in the entire assembly for every part made of this material is stored
            self.material_bounds_in_assembly: BoundBox | None = None

    def get_cms_relative_distance(
        self,
        center: tuple = (0, 0, 0)):
        """
        Calculates the relative distance between the object's center and a given center point.
        The distance is computed as the Euclidean norm between the object's center coordinates
        and the provided `center` tuple, normalized by the norm of the `center` itself.
        """
            
        try:
            relative_distance = np.linalg.norm(
                [
                    self.center.x - center[0],
                    self.center.y - center[1],
                    self.center.z - center[2],
                ]
            ) / np.linalg.norm(center)
        except ZeroDivisionError:
            relative_distance = float('inf')
        return relative_distance
        

    def get_bounding_box_relative_distance(
        self,
        bb: tuple = (0, 0, 0)):
            
        """
        Calculates the relative distance between the bounding box of the current object and a given bounding box.
        The distance is computed as the Euclidean norm of the difference between the dimensions of the current object's bounding box (`self.bb`) and the provided bounding box (`bb`), normalized by the norm of the provided bounding box dimensions.
        """
        try:
            relative_distance = np.linalg.norm(
                [self.bb.xlen - bb[0], self.bb.ylen - bb[1], self.bb.zlen - bb[2]]
            ) / np.linalg.norm(bb)
        except ZeroDivisionError:
            relative_distance = float('inf')
        return relative_distance

    def get_volume_relative_distance(self,
        volume: float = 1):
        """
        Calculate the relative distance between the object's volume and a given volume.
        Parameters:
            volume (float, optional): The reference volume to compare against. Defaults to 1.
        Returns:
            float: The absolute relative difference between the object's volume and the given volume.
        """
        
        return np.abs(self.volume - volume) / volume

    def get_similarity_score(
        self,
        center: tuple = (0, 0, 0),
        bb: tuple = (0, 0, 0),
        volume: float = 1,
        tolerance=1e-2,)->float:       
        """
        Calculates a similarity score based on the relative distances of the center of mass, bounding box, and volume 
        compared to provided reference values.
        Args:
            center (tuple, optional): Reference center of mass coordinates. Defaults to (0, 0, 0).
            bb (tuple, optional): Reference bounding box dimensions or coordinates. Defaults to (0, 0, 0).
            volume (float, optional): Reference volume. Defaults to 1.
            tolerance (float, optional): Tolerance value for similarity calculations. Defaults to 1e-2.
        Returns:
            float: The sum of the relative distances for center of mass, bounding box, and volume as a similarity score.
        """
        cms_rel_dist = self.get_cms_relative_distance(center = center)
        bb_rel_dist = self.get_bounding_box_relative_distance(bb = bb)
        vol_rel_dist = self.get_volume_relative_distance(volume = volume)
        return cms_rel_dist + bb_rel_dist + vol_rel_dist

    def object_is_similar(self, center: tuple = (0, 0, 0), bb: tuple = (0, 0, 0), volume: float = 1, tolerance=1e-2) -> bool:

        """
        Determines if the current object is similar to another based on center of mass, bounding box, and volume.

        Args:
            center (tuple, optional): The reference center of mass coordinates (x, y, z). Defaults to (0, 0, 0).
            bb (tuple, optional): The reference bounding box dimensions. Defaults to (0, 0, 0).
            volume (float, optional): The reference volume. Defaults to 1.
            tolerance (float, optional): The maximum allowed relative difference for similarity. Defaults to 1e-2.

        Returns:
            bool: True if the center of mass, bounding box, and volume are all within the specified tolerance; False otherwise.
        """
        cms_close = (
            self.get_cms_relative_distance(center = center) < tolerance
        )
        bb_close = (
            self.get_bounding_box_relative_distance(bb = bb) < tolerance
        )
        vol_close = self.get_volume_relative_distance(volume = volume) < tolerance
        return cms_close and bb_close and vol_close


    def export_stp(self):
        """
        Not yet implemented
        """
        pass


def idx_similar(entity_list, center, bounding_box, volume):
    """
    Which objects (if any) are similar to the object specified by the set of parameters

    Parameters
    ----------
    entity_list : List
        List of entities to check agains the given parameters
    center :
        The center coordinate to compare with
    bb :
        Bounding box to compare with
    volume :
        Volume to compare with

    Returns
    -------
    int
        The index in the entity_list for which a solid is similar in terms of bounding box, cms, and volume
        If no similar object is found return -1.
    """
    idx_found = []
    for i, ent in enumerate(entity_list):
        if ent.object_is_similar(
            [center.x, center.y, center.z],
            [bounding_box.xlen, bounding_box.ylen, bounding_box.zlen],
            volume,
            tolerance=1e1,
        ):
            idx_found.append(i)
    if len(idx_found) > 1:
        # we have multiple matches - pick the best one
        dsmall = 1e9
        for i, ent in enumerate([entity_list[idx] for idx in idx_found]):
            d = ent.get_similarity_score(
                [center.x, center.y, center.z],
                [bounding_box.xlen, bounding_box.ylen, bounding_box.zlen],
                volume,
            )
            if d < dsmall:
                dsmall = d
                end_idx = idx_found[i]
    elif len(idx_found) == 1:
        end_idx = idx_found[0]
        print(f"INFO: Found index at {end_idx}")
    else:
        end_idx = -1
        print("INFO: No similar object found")
    return end_idx


def similar_solids(solid1_vol, solid1_bb, solid1_c, solid2_vol, solid2_bb, solid2_c):
    """This function compares two solids and reports their similarity constant.
    defined as the sum of:
    1. cubic root difference in volume
    2. difference of bounding box diagonal
    3. difference in location vector.
    """
    dV = math.pow(
        math.fabs(solid1_vol - solid2_vol), 0.3333333333333333333333333333333333
    )
    dBB = math.fabs(solid1_bb.DiagonalLength - solid2_bb.DiagonalLength)
    c1 = solid1_c
    c2 = solid2_c
    dCntr = math.sqrt(
        (c1.x - c2.x) * (c1.x - c2.x)
        + (c1.y - c2.y) * (c1.y - c2.y)
        + (c1.z - c2.z) * (c1.z - c2.z)
    )
    return dV + dBB + dCntr


class H5MTransformer:
    """
    H5MTransformer is a class for manipulating H5M files using the pyMOAB/MOAB core library.

    Parameters:
    verbose : int
        verbosity level of output. 0: most quiet, 1: some output, 2+: a lot of diagnostic output

    Attributes:
        existing_h5m_filepath (str): Path to the original H5M file.
        updated_h5m_filepath (str): Path to the updated H5M file to be written.
        moab_core: Instance of pyMOAB/MOAB core used for file operations.

    Methods:
        __init__(h5m_filename: str):
            Initializes the transformer with the specified H5M file path.

        read_h5m_file_data() -> None:
            Loads the H5M file into the MOAB core.

        get_entity_by_id(entity_id: int):
            Retrieves an entity from the H5M file by its unique ID.

        retag_entity(material_tag: str, new_material_tag: str, entity_id: int):

        write_updated_h5m_file() -> None:
            Writes the updated H5M file to disk.

    """

    def __init__(self, h5m_filename: str, verbose: int = 1):
        self.verbose = verbose
        self.existing_h5m_filename = h5m_filename  # Path to the H5M file
        self.existing_h5m_filepath = os.path.join(os.getcwd(), h5m_filename)
        self.updated_h5m_filepath = h5m_filename.replace(".h5m", "_updated.h5m")
        self.moab_core = core.Core()  # Instance of pyMOAB / MOAB to perform operations

    def read_h5m_file_data(self) -> None:
        """
        This method reads the H5M file specified by h5m_filepath.
        """
        self.moab_core.load_file(self.existing_h5m_filepath)

    def get_entity_by_id(self, entity_id: int) -> Union[np.uint64 , None]:
        """
        This method retrieves an entity from the H5M file by its ID.
        Args:
            entity_id (int): The ID of the entity to retrieve.
        Returns:
            The entity object if found, None otherwise.
        """
        all_entities = self.moab_core.get_entities_by_handle(0)  # 0 is the root set
        for entity in all_entities:
            if entity == entity_id:
                print(type(entity))
                return entity
        return None

    def retag_entity(self, new_material_tag: str, entity_id: int):
        """
        Updates the material tag of a specified entity.
        Parameters:
            new_material_tag (str): The new material tag to assign to the entity.
            entity_id (int): The unique identifier of the entity to be retagged.
        Returns:
            None
        """
        entity = self.get_entity_by_id(entity_id)
        if entity is None:
            print(f"ERROR: Entity with ID {entity_id} not found.")
            return
        tag_handles = self.moab_core.tag_get_tags_on_entity(entity)
        for tag_handle in tag_handles:
            tag_name = tag_handle.get_name()
            if tag_name == types.NAME_TAG_NAME:
                # Update the material tag
                self.moab_core.tag_set_data(tag_handle, [entity], [[new_material_tag]])
                if self.verbose > 1:
                    print(f"INFO: Updated entity {entity_id} with new material tag '{new_material_tag}'")
                return
        if self.verbose > 0:
            print(f"WARNING: Material tag not found on entity {entity_id}. No changes made.")

    def write_updated_h5m_file(self) -> None:
        """
        This method writes the updated H5M file to disk.
        """
        self.moab_core.write_file(self.updated_h5m_filepath)

    def get_all_entity_ids(self) -> list[np.uint64]:
        """
        This method retrieves all entity IDs from the H5M file.
        Returns:
            List of entity IDs.
        """
        all_entities = self.moab_core.get_entities_by_handle(0)  # 0 is the root set
        entity_ids = list(all_entities)
        if self.verbose > 0:
            print(f"INFO: Found {len(entity_ids)} entities in the H5M file.")

        return entity_ids

    def get_all_materials(self) -> list[str]:
        """
        This method retrieves all unique material tags from the H5M file.
        Returns:
            List of unique material tags.
        """
        all_entities = self.moab_core.get_entities_by_handle(0)  # 0 is the root set
        material_tags = set()
        for entity in all_entities:
            tag_handles = self.moab_core.tag_get_tags_on_entity(entity)
            if self.verbose >= 2:
                print(
                    f"DEBUG: Entity {entity} has tags {[tag.get_name() for tag in tag_handles]}"
                )

            tag_data = {}
            for tag_handle in tag_handles:
                tag_name = tag_handle.get_name()
                if not tag_name:
                    print(
                        f"WARNING: tag_handle does not have get_name() method: {tag_handle}"
                    )
                    continue

                if tag_name == types.NAME_TAG_NAME:
                    tag_data = self.moab_core.tag_get_data(tag_handle, [entity])

                    nested_tag_list = tag_data.tolist()
                    # normalize the list of lists
                    normalized_tag_list = [
                        tag.strip()
                        for sublist in nested_tag_list
                        for tag in sublist
                        if tag.strip()
                    ]

                    for tag in normalized_tag_list:
                        material_tags.add(tag)
        if self.verbose > 0:
            print(f"INFO: Found {len(material_tags)} unique materials in the H5M file.")
        if self.verbose >= 2:
            print(f"DEBUG: Unique materials are: {material_tags}")
        return list(material_tags)

    def get_all_entities_of_material(self, material_tag: str) -> list[int]:
        """
        Retrieves all entity IDs associated with a given material tag.
        Parameters:
            material_tag (str): The material tag to search for.
        Returns:
            List of entity IDs associated with the material.
        """
        all_entities = self.moab_core.get_entities_by_handle(0)  # 0 is the root set
        matching_entity_ids = []
        for entity in all_entities:
            tag_handles = self.moab_core.tag_get_tags_on_entity(entity)
            for tag_handle in tag_handles:
                tag_name = tag_handle.get_name()
                if tag_name == types.NAME_TAG_NAME:
                    tag_data = self.moab_core.tag_get_data(tag_handle, [entity])
                    nested_tag_list = tag_data.tolist()
                    normalized_tag_list = [
                        tag.strip()
                        for sublist in nested_tag_list
                        for tag in sublist
                        if tag.strip()
                    ]
                    if material_tag in normalized_tag_list:
                        matching_entity_ids.append(entity)
                        break
        if self.verbose > 0:
            print(
                f"INFO: Found {len(matching_entity_ids)} entities for material '{material_tag}'."
            )
        return matching_entity_ids

class Assembly:
    """
    Main class representing a geometry model to process

    This class encapsulates a geometry defined by a set of step-files that are to be converted
    into a single surface-meshed geometry.
    N.b. if geometries are not overlapping it may be simple to use a single Assembly-object per
    step-file and merge them later using the merge_h5m-methiod.

    Parameters
    ----------
    stp_files : 
        List of filenames of the step-files to be processed. Most often this is a lits with single member
    verbose : int
        verbosity level of output. 0: most quiet, 1: some output, 2+: a lot of diagnostic output
    default_tag : str
        Material tag to apply to objects where none is found or where no tag-conversion is applicable
    implicit_complement : str or None
        Material tag to be applied to the volume _not_ claimed by any object. If set to None, it is ignored.

    Attributes
    ----------
    stp_files : List
        List of filenames of the step-files to be processed.
    verbose : int
        verbosity level of output. 0: most quiet, 1: some output, 2+: a lot of diagnostic output
    default_tag : string
        Material tag to apply to objects where none is found or where no tag-conversion is applicable
    delete_intermediate : bool
        Should intermediate data-files be deleted after finishing a run? Superceded by "cleanup"
    cleanup : bool
        Delete the data-directory containing the intermediate data-files?
    datadir : string
        Data directory where intermediate data is stored. If == "." (default) a unique directory will be created
    tags : dictionary
        Dictionary of material tag-conversions to be applied after extracting tags from the step-file(s). Keys
        are expected to be regexp-patterns to be matched agains the extracted tags (see also set_tag_delim), and
        values are the atcual tags to be set. E.g. tags = {"steel.*" : "niobium"} will replace all tags starting with
        "steel" with the tag "niobium".
    sequential_tags : list
        List of tages to be applied sequentially to a model. If set, this applies the tags in the list in the sequence
        they are encountered in the step-file, until either one is exhausted. If the seq. tag list is exhausted,
        remaining objects will retain the tag extracted from the step-file or get the default depending on whether
        noextract_tags is set.
    implicit_complement : string
        Material tag to be applied to the volume _not_ claimed by any object.
    noextract_tags : bool
        If true, do not extract materials tags from the step-files to fill up.
    tag_delim_pattern : string
        Regex-pattern which determines what the constitutes the first section of a part name. The first bit is the
        extracted material tag. E.g. with default setting, the names "steel pipe" and "steel_pipe" both become the tag "steel".
    """

    def __init__(
        self,
        stp_files=[],
        stl_files=[],
        verbose: int = 1,
        default_tag="vacuum",
        implicit_complement=None,
    ):
        self.stp_files = stp_files
        self.stl_files = stl_files
        self.entities = []
        self.verbose = verbose
        self.default_tag = default_tag
        self.delete_intermediate = False
        self.cleanup = False
        self.datadir="."
        self.tags = None
        self.sequential_tags = None
        self.implicit_complement = implicit_complement
        self.noextract_tags = True
        self.tag_delim_pattern=r"^([^\s_@]+)"

        #check if we can write a moab-detabase.
        self.dummy_h5m()

    def dummy_h5m(self):
        mbc, mbt = self.init_moab()
        mbc.write_file("dummy.h5m")
        try:
            self.check_h5m_file("dummy.h5m")
        except RuntimeError as e:
            raise e
        os.unlink("dummy.h5m")
        return True

    def set_tag_delim(self,delimiters: str):
        self.tag_delim_pattern=r"^([^" + delimiters + r"]+)"

    def run(
        self,
        backend: str = "stl",
        h5m_filename: str = "dagmc.h5m",
        merge: bool = True,
        imprint: bool = False,
        **kwargs: dict,
    ):
        """
        Convenience function to easily run the processing sequence

        This method merges the calls to import_stp_files, imprint, merge, and solids_to_h5m into a single call

        parameters
        ----------
        backend : string
            Mesh creatinon backend to use. Allowed values are 'stl', 'stl2', 'gmsh, and 'db'. The currently 
            recommended choices ae "stl2" and "db".
        h5m_filename : string
            Filename of the created output data-file which will contains the meshed geometry
        merge: bool
            Perform a merge step during geometry processing
        imprint: bool
            Perform imprinting step during geometry processing

        """
        self.import_stp_files(**kwargs)
        if imprint:
            self.imprint_all()
        if merge:
            self.merge_all()
        self.solids_to_h5m(backend=backend, h5m_filename=h5m_filename, **kwargs)

    def import_stp_files(
        self,
        tags: dict = None,
        sequential_tags: iter = None,
        match_anywhere: bool = False,
        default_tag: str = "vacuum",
        scale: float = 0.1,
        translate: iter = [],
        rotate: iter = [],
        vol_skip: iter = [],
        **kwargs: dict,
    ):
        """
        Import a list of step-files and extract material tags.

        parameters
        ----------
        tags: dictionary
            Contains pairs of reg.exp. patterns and material tags. If not None, entities with
            names matching the patterns will be assigned to corresponding tag. If no patterns match
            the default_tag will be applied
        match_anywhere: bool
            Match patterns anywhere in the entitiy name
        default_tag: string
            The material tag that will be applied if no patterns are matched.
        scale: float
            Overall scaling factor applied to all parts
        translate: List
            3D translation vector to apply to all parts in the step-file.
        rotate: list
            Rotation angles to apply to the parts in the step-file.
        vol_skip: list
            Numbers of volumes to skip meshing for.
        """
        for stp in self.stp_files:
          warn, ct = has_degenerate_toroids(stp,True)
          if warn:
            print(f'WARNING: Step file {stp} has {ct} degenerate toroid surfaces. These are known to cause problems in some cases',file=sys.stderr)

        tags_set = 0
        # clear list to avoid double-import
        self.entities = []

        message = "Need gmsh python module installed to extract material tags from step-file. please supply a 'sequential_tags'-list instead"
        assert (nogmsh and tags is None and sequential_tags is not None) or (
            not nogmsh
        ), message
        # if no gmsh module was imported we must rely on explicit sequential tags, so check they're there.
        i = 1
        for stp in self.stp_files:
            solid = self.load_stp_file(stp, scale, translate, rotate)

            ents = []
            # try if solid is iterable
            try:
                for j in range(len(solid)):
                    if j + i not in vol_skip:
                        s = solid[j]
                        e = Entity(solid=s)
                        ents.append(e)
                i = i + len(solid)
            except Exception as _e:
                e = Entity(solid=solid)
                ents.append(e)

            if tags is None and sequential_tags is None:
                # also import using gmsh to extract the material tags from the labels in the step files
                gmsh.initialize()
                vols = gmsh.model.occ.importShapes(stp)
                gmsh.model.occ.synchronize()
                for e, v in zip(ents, vols):
                    vid = v[1]
                    try:
                        s = gmsh.model.getEntityName(3, vid)
                        part = s.split("/")[-1]
                        g = re.match(self.tag_delim_pattern, part)
                        tag = unidecode(g[0])
                        if self.verbose > 1:
                            print(
                                f"INFO: Tagging volume #{vid} label:{s} with material {tag}"
                            )
                    except Exception as _e:
                        tag = default_tag
                    e.tag = tag
                    tags_set = tags_set + 1
                gmsh.finalize()
            elif tags:
                # tag objects according to the tags dictionary.
                gmsh.initialize()
                vols = gmsh.model.occ.importShapes(stp)
                gmsh.model.occ.synchronize()
                for e, v in zip(ents, vols):
                    vid = v[1]
                    s = gmsh.model.getEntityName(3, vid)
                    part = s.split("/")[-1]
                    tag = None
                    for k in tags.keys():
                        if match_anywhere:
                            g = re.search(k, part)
                        else:
                            g = re.match(k, part)
                        if g is not None:
                            tag = tags[k]
                            tags_set = tags_set + 1
                            break
                    #if tag is still not set at this point we will either leave it or set it to the default.
                    if tag is None:
                        if self.noextract_tags:
                            tag = self.default_tag
                        else:
                            #use tag from stepfile
                            try:
                                g = re.match(self.tag_delim_pattern, part)
                                tags_set = tags_set + 1
                                tag=g[0]
                            except:
                                #this e.g. happens when there is no tag in the step-file
                                tag=self.default_tag
                    #apply the selected tag to the entity
                    e.tag = tag
                    if self.verbose > 1:
                        print(
                            f"INFO: Tagging volume #{vid} label:{s} with material {tag}"
                        )
                gmsh.finalize()
            elif sequential_tags:
                for e, t in zip_longest(
                    ents, sequential_tags[tags_set:], fillvalue=self.default_tag
                ):
                    # Apply tags to the imported volumes in the sequence they get imported.
                    if e == default_tag:
                        # this means we have exhausted the ents list
                        break
                    e.tag = t
                tags_set += len(ents)

            self.entities.extend(ents)
        if tags_set != len(self.entities):
            print(
                f"WARNING: {len(self.entities)-tags_set} volumes were tagged with the default ({default_tag}) material."
            )
        # Iterate over the ents array and get all unique tags
        unique_tags = set(e.tag for e in self.entities if e.tag != default_tag)
        if unique_tags:
            print(f"INFO: Found unique tags: {unique_tags}")
        # Iterate over unique tags
        for tag in unique_tags:
            print(f"INFO: Processing tag: {tag}")
            # Get all entities with this tag
            tagged_entities = [e for e in self.entities if e.tag == tag]
            print(f"INFO: Found {len(tagged_entities)} entities with tag {tag}")
            # Get the bounding box for this tag
            if tagged_entities:
                # get the solid from the tagged_entities

                solid_bounding_box_array = [
                    e.solid.BoundingBox() for e in tagged_entities
                ]
                all_solids_bounding_box = solid_bounding_box_array[0]
                for bb in solid_bounding_box_array[1:]:
                    all_solids_bounding_box = all_solids_bounding_box.add(bb)
                bbox = all_solids_bounding_box
                bbox_text = f"xmin={bbox.xmin}, xmax={bbox.xmax}, ymin={bbox.ymin}, ymax={bbox.ymax}, zmin={bbox.zmin}, zmax={bbox.zmax}"
                print(f"INFO: Bounding box for tag {tag}: {bbox_text}")
                # Assign this bounding box to the material_bounds_in_assembly attribute
                for e in tagged_entities:
                    e.material_bounds_in_assembly = bbox

    def load_stp_file(
        self,
        filename: str,
        scale_factor: float = 0.1,
        translate: list = [],
        rotate: list = [],
    ):
        """Loads a stp file and makes the 3D solid and wires available for use.

        Args:
            filename: the filename of the file containing the step-definition.
            scale_factor: a scaling factor to apply to the geometry that can be
                used to increase the size or decrease the size of the geometry.
                Useful when converting the geometry to cm for use in neutronics
                simulations. The default (0.1) corresponds to the standard setting of OpenCASCADE
                which assumes mm to the the standard of OpenMC which is cm.
            translate: optional, iterable of 3 float iterable. Translation vector to
                apply to parts in the model before meshing occurs. If list of list then
                each translation in the list is applied to 1 part in sequence.

        Returns:
            [CadQuery.solid]
        """
        # import _all_ the shapes in the file - i.e. may return a list
        part = cq.importers.importStep(str(filename)).vals()

        # apply apply_transforms
        scaled_part = self.apply_transforms(
            part, filename, scale_factor, translate, rotate
        )

        solid = []
        # serialize
        # Solids() returns a list even if the part is not a Compund object
        try:
            for p in scaled_part:
                solid.extend(p.Solids())
        except Exception as _e:
            solid.extend(scaled_part.Solids())

        return solid

    def apply_transforms(self, part, filename, scale_factor, translate, rotate):
        # scale the shapes even if the factor is 1.
        if self.verbose != 0:
            print(f"INFO: {str(filename)} imported - scaling")
        try:
            transformed_part = [p.scale(scale_factor) for p in part]
        except Exception as _e:
            transformed_part = part.scale(scale_factor)

        # translation
        if translate:
            if self.verbose != 0:
                print(f"INFO: {str(filename)} imported - applying translation(s)")
            try:
                vols = translate[::2]
                translations = translate[1::2]
                for v in enumerate(vols):
                    if self.verbose > 1:
                        print(
                            f"INFO: Applying translation: {translations[v[0]]} to vol(s) {v[1]}"
                        )
                    for vol in v[1]:
                        transformed_part[vol - 1] = transformed_part[vol - 1].translate(
                            translations[v[0]]
                        )
            except Exception as _e:
                transformed_part = transformed_part.translate(translate[1])

        # rotation
        if rotate:
            if self.verbose != 0:
                print(f"INFO: {str(filename)} imported - applying rotation(s)")
            try:
                vols = rotate[::4]
                for v in enumerate(vols):
                    if self.verbose > 1:
                        print(
                            f"INFO: Applying rotation: {rotate[v[0]+3]} degrees about ax {rotate[v[0]+1]},{rotate[v[0]+2]} to vol(s) {v[1]}\n"
                        )
                    for vol in v[1]:
                        transformed_part[vol - 1] = transformed_part[vol - 1].rotate(
                            rotate[4 * v[0] + 1],
                            rotate[4 * v[0] + 2],
                            rotate[4 * v[0] + 3],
                        )
            except Exception as _e:
                transformed_part = transformed_part.rotate(
                    rotate[1], rotate[2], rotate[3]
                )

        return transformed_part

    def print_summary(self):
        # output a summary of the meshing results
        print(f'SUMMARY: {"solid_id":8} {"material_tag":16} {"stl-file":16}')
        for i, a in zip(range(len(self.entities)), self.entities):
            if not isinstance(a.stl, str):
                print(
                    f"SUMMARY: {i+1:8} {a.tag:16} "
                    + " ".join([f"{stl[0]:16}" for stl in a.stl])
                )
            else:
                print(f"SUMMARY: {i+1:8} {a.tag:16} {a.stl:16}")

    def _datadir_name(self,h5m_filename=""):
        h5mf=pl.Path(h5m_filename)
        if self.datadir==".":
            datadir = datetime.now().strftime(f"{h5mf.stem}_%Y%m%d_%H%M%S.%f")
        else:
            datadir = self.datadir
        if (self.verbose):
            print(f"INFO: storing temporary data in directory: {datadir}")
        return datadir

    def solids_to_h5m(
        self,
        brep_filename: str = None,
        h5m_filename: str = "dagmc.h5m",
        samples: int = 100,
        backend: str = "gmsh",
        heal: bool = True,

        **kwargs: dict,
    ):
        """
        Convert the imported list of solids into surface-meshed bounded objects

        Main processing function that processes the list of already imported entities. It creates a
        mesher backend object as appropriate and calls that to create a surface mesh.

        parameters
        ----------
        h5m_filename: string
            Name of the output file
        backend: string
            Name of the meshing backend to call. Allowed values are in ["gmsh", "stl", "stl2", "db"]
        heal: bool
            Use trimesh (if available) to perform a healing step on the surface normals"
        """
        #if h5m_filename is not in cwd, run where it resides.
        cwd=pl.Path.cwd()
        h5m_path=pl.Path(h5m_filename)
        os.chdir(h5m_path.parent)
        with mesher_datadir(self._datadir_name(h5m_path.name),self.cleanup, movein=True) as datadir:
            mesher_config["entities"] = self.entities
            meshgen = am.meshers.get(backend, **mesher_config)
            meshgen.set_verbosity(self.verbose)
            stl_list = meshgen.generate_stls()

            for e, s in zip(self.entities, stl_list):
                e.stl = s
            if self.verbose:
                self.print_summary()
            if backend in [ "stl2", "db" ] :
                self.stl2h5m_byface(h5m_path.name, True)
            else:
                if heal:
                    stl_list = self.heal_stls(stl_list)
                self.stl2h5m(stl_list, h5m_path.name, True)
        os.chdir(cwd)

    def stl2h5m_byface(self, h5m_file: str = "dagmc.h5m", vtk: bool = False) -> str:
        """create a h5m-file with a moab structure and fills
        it with the dagmc structure using the pymoab framework. Optionally creates
        a vtk-file.
        """
        if self.verbose > 0:
            print("INFO: reassembling stl-files into h5m structure")
        h5m_p = pl.Path(h5m_file)
        mbcore, mbtags = self.init_moab()
        mbcore = self.add_entities_to_moab_core(mbcore, mbtags)

        all_sets = mbcore.get_entities_by_handle(0)
        file_set = mbcore.create_meshset()

        mbcore.add_entities(file_set, all_sets)

        if self.verbose > 0:
            print(f'INFO: writing geometry to h5m: "{h5m_file}".')
        mbcore.write_file(str(h5m_p))

        self.check_h5m_file(h5m_file)
        if vtk:
            if self.verbose > 0:
                print(f'INFO: writing geometry to vtk: ' + str(h5m_p.with_suffix(".vtk")) )
            mbcore.write_file(str(h5m_p.with_suffix(".vtk")))

        return str(h5m_p)

    def stl2h5m(
        self, stls: list, h5m_file: str = "dagmc.h5m", vtk: bool = False
    ) -> str:
        """function that exports the list of stls that we have presumably generated somehow
        and merges them into a DAGMC h5m-file by means of the MOAB-framework.
        """

        if self.verbose > 0:
            print("INFO: reassembling stl-files into h5m structure")
        h5m_p = pl.Path(h5m_file)
        moab_core, moab_tags = self.init_moab()

        sid, vid = (1, 1)
        is_last=False
        for e in self.entities:
            if self.verbose > 1:
                print(
                    f'INFO: add stl-file "{e.stl}" with tag "{e.tag}" to MOAB structure'
                )
            if (e == self.entities[-1]):
                is_last=True

            moab_core = self.add_stl_to_moab_core(
                moab_core, sid, vid, e.tag, moab_tags, e.stl,
                last=is_last
            )
            vid += 1
            sid += 1
            if self.delete_intermediate:
                p = pl.Path(e.stl)
                p.unlink()

        all_sets = moab_core.get_entities_by_handle(0)

        file_set = moab_core.create_meshset()

        moab_core.add_entities(file_set, all_sets)

        if self.verbose > 0:
            print(f'INFO: writing geometry to h5m: "{h5m_file}".')
        moab_core.write_file(str(h5m_p))

        self.check_h5m_file(h5m_file)

        if vtk:
            moab_core.write_file(str(h5m_p.with_suffix(".vtk")))

        return str(h5m_p)

    def check_h5m_file(self, h5m_file: str = "dagmc.h5m"):
        with open(h5m_file, "rb") as f:
            magic_bytes = f.read(8)
            if magic_bytes != b"\x89HDF\x0d\x0a\x1a\x0a":
                raise RuntimeError(
                    "Generated h5m-file does not appear to be a hdf-file"
                )

    def add_entities_to_moab_core(self, mbcore: core.Core, mbtags: dict, noimplicit=False, in_datadir="."):
        vsets = []
        glob_id = 0
        for i in range(len(self.entities)):
            vset = mbcore.create_meshset()
            vsets.append(vset)
            glob_id += 1
            mbcore.tag_set_data(mbtags["global_id"], vset, glob_id)
            mbcore.tag_set_data(mbtags["geom_dimension"], vset, 3)
            mbcore.tag_set_data(mbtags["category"], vset, "Volume")

        faces_added = {}
        sid = 0
        gid = 0
        for i, e in enumerate(self.entities):
            for j, T in enumerate(e.stl):
                f, sense = T
                if f not in faces_added:
                    fset = mbcore.create_meshset()
                    sid += 1
                    glob_id += 1
                    faces_added[f] = fset

                    mbcore.tag_set_data(mbtags["global_id"], fset, glob_id)
                    mbcore.tag_set_data(mbtags["geom_dimension"], fset, 2)
                    mbcore.tag_set_data(mbtags["category"], fset, "Surface")

                    mbcore.add_parent_child(vsets[i], fset)
                    if len(sense) == 2 and sense[1]!=-1:
                        mbcore.tag_set_data(
                            mbtags["surf_sense"],
                            fset,
                            np.array(
                                [vsets[sense[0]], vsets[sense[1]]], dtype="uint64"
                            ),
                        )
                    else:
                        mbcore.tag_set_data(
                            mbtags["surf_sense"],
                            fset,
                            np.array([vsets[sense[0]], 0], dtype="uint64"),
                        )
                    mbcore.load_file(str(pl.Path(in_datadir) / f), fset)
                else:
                    # this face has already been added so only add a parent child relation here
                    fset = faces_added[f]
                    mbcore.add_parent_child(vsets[i], fset)
            # make this a group, this could ideally be a set of volumes with the same material
            gset = mbcore.create_meshset()
            gid += 1
            glob_id += 1
            mbcore.tag_set_data(mbtags["category"], gset, "Group")
            # reflective is a special case that should not have mat: in front
            if not e.tag == "reflective":
                dagmc_material_tag = f"mat:{e.tag}"
            else:
                dagmc_material_tag = e.tag

            mbcore.tag_set_data(mbtags["name"], gset, dagmc_material_tag)
            mbcore.tag_set_data(mbtags["geom_dimension"], gset, 4)

            # add the volume to this group set
            mbcore.add_entity(gset, vsets[i])

        # if wanted add an implicit complement material
        if ( self.implicit_complement is not None and not noimplicit):
            gset = mbcore.create_meshset()
            mbcore.tag_set_data(mbtags["category"], gset, "Group")
            mbcore.tag_set_data(mbtags["name"], gset, f"mat:{self.implicit_complement}_comp")
            mbcore.tag_set_data(mbtags["geom_dimension"], gset, 4)
            mbcore.add_entity(gset, vsets[-1])

        # finally set the faceting tolerance tag
        mbcore.tag_set_data(
            mbtags["faceting_tol"], mbcore.get_root_set(), np.array((mesher_config["tolerance"],))
        )

        return mbcore

    def add_stl_to_moab_core(
        self,
        mbcore: core.Core,
        surface_id: int,
        volume_id: int,
        material_name: str,
        mbtags: dict,
        stl_filename: str,
        last: bool = False
    ) -> core.Core:
        """
        Appends a set of surfaces (comprising a volume) from an stl-file to a moab.Core object and returns the updated object

        Args:
            mbcore: A moab core object
            surface_id: the id number to apply to the surface
            volume_id: the id numbers to apply to the volumes
            material_name: the material tag name to add. the value provided will
                be prepended with "mat:" unless it is "reflective" which is
                a special case and therefore will remain as is.
            mbtags: A dictionary of the MOAB tags
            stl_filename: the filename of the stl file to load into the moab core

        Returns:
            An updated pymoab.core.Core() instance
        """

        surface_set = mbcore.create_meshset()
        volume_set = mbcore.create_meshset()

        # recent versions of MOAB handle this automatically
        # but best to go ahead and do it manually
        mbcore.tag_set_data(mbtags["global_id"], volume_set, volume_id)
        mbcore.tag_set_data(mbtags["global_id"], surface_set, surface_id)

        # set geom IDs
        mbcore.tag_set_data(mbtags["geom_dimension"], volume_set, 3)
        mbcore.tag_set_data(mbtags["geom_dimension"], surface_set, 2)
        # set category tag values
        mbcore.tag_set_data(mbtags["category"], volume_set, "Volume")
        mbcore.tag_set_data(mbtags["category"], surface_set, "Surface")

        # establish parent-child relationship
        mbcore.add_parent_child(volume_set, surface_set)

        # Set surface sense
        # This should be fixed - we should know which volume comes next, instead of just setting it to be 0
        sense_data = [volume_set, np.uint64(0)]
        mbcore.tag_set_data(mbtags["surf_sense"], surface_set, sense_data)

        # load the stl triangles/vertices into the surface set
        mbcore.load_file(stl_filename, surface_set)

        group_set = mbcore.create_meshset()
        mbcore.tag_set_data(mbtags["category"], group_set, "Group")

        # reflective is a special case that should not have mat: in front
        if not material_name == "reflective":
            dag_material_tag = "mat:{}".format(material_name)
        else:
            dag_material_tag = material_name

        mbcore.tag_set_data(mbtags["name"], group_set, dag_material_tag)

        mbcore.tag_set_data(mbtags["geom_dimension"], group_set, 4)

        # add the volume to this group set
        mbcore.add_entity(group_set, volume_set)

        # if this is the last element also add some extra information
        if ( last ):
            # if wanted add an implicit complement material
            if ( self.implicit_complement is not None ):
                gset = mbcore.create_meshset()
                mbcore.tag_set_data(mbtags["category"], gset, "Group")
                mbcore.tag_set_data(mbtags["name"], gset, f"mat:{self.implicit_complement}_comp")
                mbcore.tag_set_data(mbtags["geom_dimension"], gset, 4)
                mbcore.add_entity(gset, volume_set)

            # finally set the faceting tolerance tag
            mbcore.tag_set_data(
                mbtags["faceting_tol"], mbcore.get_root_set(), np.array((mesher_config["tolerance"],))
            )

        return mbcore

    def init_moab(self):
        """Creates a MOAB Core instance which can be built up by adding sets of
        triangles to the instance
        Returns:
        (pymoab Core): A pymoab.core.Core() instance
        (pymoab tag_handle): A pymoab.core.tag_get_handle() instance
        """
        # create pymoab instance
        moab_core = core.Core()

        tags = dict()
        sense_tag_name = "GEOM_SENSE_2"
        sense_tag_size = 2
        tags["surf_sense"] = moab_core.tag_get_handle(
            sense_tag_name,
            sense_tag_size,
            types.MB_TYPE_HANDLE,
            types.MB_TAG_SPARSE,
            create_if_missing=True,
        )
        tags["category"] = moab_core.tag_get_handle(
            types.CATEGORY_TAG_NAME,
            types.CATEGORY_TAG_SIZE,
            types.MB_TYPE_OPAQUE,
            types.MB_TAG_SPARSE,
            create_if_missing=True,
        )
        tags["name"] = moab_core.tag_get_handle(
            types.NAME_TAG_NAME,
            types.NAME_TAG_SIZE,
            types.MB_TYPE_OPAQUE,
            types.MB_TAG_SPARSE,
            create_if_missing=True,
        )
        tags["geom_dimension"] = moab_core.tag_get_handle(
            types.GEOM_DIMENSION_TAG_NAME,
            1,
            types.MB_TYPE_INTEGER,
            types.MB_TAG_DENSE,
            create_if_missing=True,
        )
        tags["faceting_tol"] = moab_core.tag_get_handle(
            "FACETING_TOL",
            1,
            types.MB_TYPE_DOUBLE,
            types.MB_TAG_DENSE,
            create_if_missing=True,
        )
        # Global ID is a default tag, just need the name to retrieve
        tags["global_id"] = moab_core.tag_get_handle(types.GLOBAL_ID_TAG_NAME)
        return moab_core, tags

    def merge_all(self):
        # merging a single object does not really make sense
        if len(self.entities) > 0:
            # extract cq solids backend algorithm
            unmerged = [e.solid for e in self.entities]

            # Pre-calculate and cache the volume, bounding box, and center of each
            unmerged_vols = [x.Volume() for x in unmerged]
            unmerged_bbs = [x.BoundingBox() for x in unmerged]
            unmerged_centers = [x.Center() for x in unmerged]

            # do merge
            # merged = self._merge_solids(unmerged, fuzzy_value=1e-6)
            merged = self._merge_solids_full(unmerged, fuzzy_value=1e-6)
            # the merging process may result in extra volumes.
            # We need to make sure these are at the end of the list
            # If not this results in a loss of volumes in the end.
            print("INFO: reordering volumes after merge")
            tmp_ents = []

            # Figure of which of the merged solids best corresponds to
            # each of the unmerged volumes.
            merged_solids = merged.Solids()

            # Pre-calculate and cache the volume, bounding box, and center of each
            merged_vols = [x.Volume() for x in merged_solids]
            merged_bbs = [x.BoundingBox() for x in merged_solids]
            merged_centers = [x.Center() for x in merged_solids]

            for j, orig in enumerate(unmerged):
                d_small = 1e9
                i_small = -1
                if self.verbose > 1:
                    print(
                        f"INFO: {len(merged_solids)} merged solids left in list of originally {len(merged.Solids())}"
                    )
                for i, ms in enumerate(merged_solids):
                    d = similar_solids(
                        unmerged_vols[j],
                        unmerged_bbs[j],
                        unmerged_centers[j],
                        merged_vols[i],
                        merged_bbs[i],
                        merged_centers[i],
                    )
                    if d < d_small:
                        i_small, d_small = i, d
                if i_small == -1:
                    print(
                        f"WARNING: Could not find a matching merged volume for volume {j+1}.",
                        end=" ",
                    )
                    print(
                        "This volume/entity will be skipped. Please examine the output volume carefully."
                    )
                else:
                    # Transfer the picked solid to the list of merged solids, removing (pop) it from the list
                    # This to avoid going through the whole list more than once.
                    ent = self.entities[j]
                    ent.solid = merged_solids.pop(i_small)
                    merged_vols.pop(i_small)
                    merged_bbs.pop(i_small)
                    merged_centers.pop(i_small)
                    tmp_ents.append(ent)
            self.entities = tmp_ents

    def imprint_all(self):
        if len(self.entities) > 0:
            unmerged = [e.solid for e in self.entities]
            merged = self.imprint_solids(unmerged)
            for e, m in zip(self.entities, merged):
                e.solid = m

    def _merge_solids_full(self, solids, fuzzy_value=1e-6):
        """
        adds all solids on the solids-vector to the algorithm and
        merges them as arguments. This is stabler, but more resource
        intensive
        """
        bldr = OCP.BOPAlgo.BOPAlgo_MakeConnected()
        for shape in solids:
            if isinstance(shape, cq.occ_impl.shapes.Compound):
                bldr.AddArgument(shape.wrapped)
            else:
                try:
                    bldr.AddArgument(shape.val().wrapped)
                except Exception as _e:
                    bldr.AddArgument(shape.wrapped)

        bldr.SetRunParallel(False)
        bldr.Perform()

        if self.verbose > 1:
            print("INFO: Generate compound shape")
        merged = cq.Compound(bldr.Shape())
        return merged

    def _merge_solids(self, solids, fuzzy_value):
        """merge a set of cq-solids
        returns as cq-compound object
        """
        bldr = OCP.BOPAlgo.BOPAlgo_Splitter()
        bldr.SetFuzzyValue(fuzzy_value)
        # loop through all objects in geometry and split and merge them accordingly
        # shapes should be a compund cq object or a list thereof
        for i, shape in enumerate(solids):
            if self.verbose:
                print(f"splitting obj {i} of {len(solids)}")
            # checks if solid is a compound as .val() is not needed for compunds
            if isinstance(shape, cq.occ_impl.shapes.Compound):
                bldr.AddArgument(shape.wrapped)
            else:
                try:
                    bldr.AddArgument(shape.val().wrapped)
                except Exception as _e:
                    bldr.AddArgument(shape.wrapped)
            for j, shape2 in enumerate(solids):
                if i == j:
                    # don't split an object with itself
                    continue
                print(f"splitting object {i} with object {j}")
                if isinstance(shape2, cq.occ_impl.shapes.Compound):
                    bldr.AddArgument(shape2.wrapped)
                else:
                    try:
                        bldr.AddArgument(shape2.val().wrapped)
                    except Exception as _e:
                        bldr.AddArgument(shape2.wrapped)
                bldr.Perform()

        if self.verbose > 1:
            print("INFO: Commence perform step of merge")
        bldr.Perform()

        if self.verbose > 1:
            print("INFO: Commence image step of merge")
        bldr.Images()

        if self.verbose > 1:
            print("INFO: Generate compound shape")
        merged = cq.Compound(bldr.Shape())

        return merged

    def imprint_solids(self, solids):
        merged = solids
        for i in range(len(merged)):
            s0 = merged[i]
            if self.verbose > 0:
                print(f"INFO: Self_imprinting {i}")
            s0 = self.imprint_solid_on_self(s0)
            merged[i] = s0

        for i in range(len(merged)):
            if self.verbose > 0:
                print(f"INFO: Imprint on solid {i}")
            s0 = merged[i]
            for j in range(0, len(merged)):
                if j == i:
                    continue
                    # imprinting solid on itself - meaning we imprint its faces on each other.
                    s0 = self.imprint_solid_on_self(s0)
                    s0.mesh(mesher_config["tolerance"])
                    merged[i] = s0
                else:
                    s1 = merged[j]
                if self.verbose > 1:
                    print(f"INFO: Imprinting solid {j} on {i}")
                result = self.imprint_solid_on_solid(s0, s1)
                if result is None:
                    if self.verbose > 1:
                        print(
                            f"INFO: solid {j}'s bounding box does not touch/overlap that of solid {i}. Skipping imprinting"
                        )
                    merged[i] = s0
                    merged[j] = s1
                else:
                    if self.verbose > 1:
                        print(
                            f"INFO: solid {j} is possibly connected to solid {i} - perform imprint."
                        )
                    s0_i = result
                    for k, idx in enumerate([i, j]):
                        try:
                            merged[idx] = s0_i.Solids()[k]
                        except Exception as _e:
                            merged[idx] = merged[idx]
        return merged

    def imprint_solid_on_self(self, s0):
        bldr = OCP.BOPAlgo.BOPAlgo_MakeConnected()
        for fc in s0.Faces():
            bldr.AddArgument(fc.wrapped)
        bldr.SetUseOBB(True)
        bldr.Perform()
        s1 = cq.Shape(bldr.Shape())
        return s1

    def overlap_bounding_boxes(self, solid1, solid2):
        (bb1, bb2) = (solid1.BoundingBox(), solid2.BoundingBox())
        outside = (
            (bb2.xmin > bb1.xmax)
            or (bb2.xmax < bb1.xmin)
            or (bb2.ymin > bb1.ymax)
            or (bb2.ymax < bb1.ymin)
            or (bb2.zmin > bb1.zmax)
            or (bb2.zmax < bb1.zmin)
        )
        return not outside

    def imprint_solid_on_solid(self, solid0, solid1):
        # if the bounding boxes of the solids don't overlap
        # don't do anything
        if not self.overlap_bounding_boxes(solid0, solid1):
            return None

        bldr = OCP.BOPAlgo.BOPAlgo_MakeConnected()
        if isinstance(solid0, cq.occ_impl.shapes.Compound):
            bldr.AddArgument(solid0.wrapped)
        else:
            try:
                bldr.AddArgument(solid0.val().wrapped)
            except Exception as _e:
                bldr.AddArgument(solid0.wrapped)
        if isinstance(solid1, cq.occ_impl.shapes.Compound):
            bldr.AddArgument(solid1.wrapped)
        else:
            try:
                bldr.AddArgument(solid1.val().wrapped)
            except Exception as _e:
                bldr.AddArgument(solid1.wrapped)
        bldr.SetUseOBB(True)
        bldr.Perform()
        return cq.Compound(bldr.Shape())

    def heal_stls(self, stls):
        #simply return early if trimesh is not available
        try:
            import trimesh
        except ImportError:
            return

        if self.verbose > 0:
            print("INFO: checking surfaces and reparing normals")
        for e in self.entities:
            stl = e.stl
            mesh = trimesh.load_mesh(stl)
            if self.verbose > 1:
                print("INFO: stl-file", stl, ": mesh is watertight", mesh.is_watertight)
            trimesh.repair.fix_normals(
                mesh
            )  # reqired as gmsh stl export from brep can get the inside outside mixed up
            new_filename = stl[:-4] + "_with_corrected_face_normals.stl"
            mesh.export(new_filename)
            if self.delete_intermediate:
                p = pl.Path(e.stl)
                p.unlink()
            e.stl = new_filename

    def get_all_tags(self):
        # extract_all_tags from the list of entities
        return [e.tag for e in self.entitites]

    def get_unique_tags(self):
        # extract a set of unique tags
        return {self.get_all_tags()}

def merge2h5m(assemblies =[], h5m_file: str ="dagmc.h5m", vtk: bool = True, verbose: int = 1):
    """
    Function that (re)performs the assembly of an h5m_file from a set of already triangularized
    assemblies.

    Usage:
        a=Assembly(['stepfileA.step'])
        a.run( ... )
        b=Assembly(['stepfileB.step'])
        b.run( ... )
        merge2h5m([a,b],h5m_file='c.h5m')

    paramters
    ---------
    assemblies: list
        Iterable of Assemblies.
    h5m_file: string
        Filename of merged file.
    vtk: bool
        If True, also write a vtk-file of the structure.
    verbose: int
        If == 0 do not write status messages to console.
    """

    #create a dummy object - this will not actually be used for anything.
    amb=assemblies[0]
    if verbose > 0:
        print(f"INFO: reassembling stl-files into h5m structure {h5m_file}")
    h5m_p = pl.Path(h5m_file)
    mbcore, mbtags = amb.init_moab()
    all_sets = mbcore.get_entities_by_handle(0)

    for a in assemblies[:-1]:
        mbcore = a.add_entities_to_moab_core(mbcore, mbtags, noimplicit=True, in_datadir=a.datadir)
    assemblies[-1].add_entities_to_moab_core(mbcore, mbtags, noimplicit=False, in_datadir=assemblies[-1].datadir)

    if verbose > 0:
        print(f'INFO: writing geometry to h5m: "{h5m_file}".')
    mbcore.write_file(str(h5m_p))

    amb.check_h5m_file(h5m_file)
    if vtk:
        if verbose > 0:
            print(f'INFO: writing geometry to vtk: ' + str(h5m_p.with_suffix(".vtk")) )
        mbcore.write_file(str(h5m_p.with_suffix(".vtk")))
    
def hdf5_in_moab(cls):
    """
    function to perform a check on whether the written h5m-name file is actually hdf5

    Usage:
        a=CAD_to_OpenMC.assembly.hdf5_in_moab()
    """
    a = Assembly()
    try:
        a.dummy_h5m()
    except RuntimeError as e:
        print(
            "Warning: Can't write a hdf-file. Did you compile moab with hdf5 support?"
            "The resulting meshed file will actually be a vtk-file."
            f"{e}"
        )
    return True

