import sys
import logging
import os
from os.path import basename,join,splitext,exists
from collections import defaultdict


prefix = "prefix"

# You probably don't want to change anything below
working_dir = join(os.getcwd() + os.sep, prefix) # Current directory + prefix folder
mesh_folder = "mesh"
output_folder = "out"
default_bc = "18 Q 0"
RWK = 300000000
IWK = 34000000


class ExtendedSlice():
    '''Object accept extended slice notation (for example: ExtendedSlice()[:3,3:6,6:]) and
       return generator object that generates slice() objects'''
    def __getitem__(self, extendedSlice):
        for _slice in extendedSlice:
            yield _slice


class PgstkConfig():
    def __init__(self, prefix, working_dir, mesh_folder, output_folder):
        self.prefix = prefix
        self.working_dir = working_dir
        self.mesh_dir = join(working_dir + os.sep, mesh_folder)
        self.out_dir = join(working_dir + os.sep, output_folder)
    
    def prepare_directories(self):
        if not exists(self.mesh_dir):
            os.makedirs(self.mesh_dir)
        if not exists(self.out_dir):
            os.makedirs(self.out_dir)
    
    def set_blk_precision(self, nblocks):
        self.blk_precision = len(str(nblocks))
    
    def set_default_bc(self, default_bc):
        self.default_bc = default_bc
    
    def set_memory_vectors(self, RWK, IWK):
        self.RWK = RWK
        self.IWK = IWK


def get_multiple_slices(listInstance, slices):
    '''Argument "slices" must be generator object that returns slice() objects'''
    for _slice in slices:
        yield listInstance[_slice]


def cells_to_vrt(cells_range, patch_dir):
    '''cells_range is patch range from cfx4 file;
       patch_dir is the direction, from the cfx manual
       0 = solid (3-D patch),
       1 = high i, 2 = high j, 3 = high k
       4 = low i, 5 = low j, 6 = low k'''
    vrt_range = cells_range[:] # We won't change initial range

    # high i
    if patch_dir == 1:
        vrt_range[0] += 1
        vrt_range[1] += 1
        vrt_range[3] += 1
        vrt_range[5] += 1
    # high j
    elif patch_dir == 2:
        vrt_range[1] += 1
        vrt_range[2] += 1
        vrt_range[3] += 1
        vrt_range[5] += 1
    # high k
    elif patch_dir == 3:
        vrt_range[1] += 1
        vrt_range[3] += 1
        vrt_range[4] += 1
        vrt_range[5] += 1
    # low i
    elif patch_dir == 4:
        vrt_range[3] += 1
        vrt_range[5] += 1
    # low j
    elif patch_dir == 5:
        vrt_range[1] += 1
        vrt_range[5] += 1
    # low k
    elif patch_dir == 6:
        vrt_range[1] += 1
        vrt_range[3] += 1
    
    return vrt_range


def invert_indexes(array, index1, index2):
    array[index1], array[index2] = array[index2], array[index1]
    

def correct_patch_range(patch_range, dirIndex1, dirIndex2, dirIndex3):
    '''patch_range is patch range from cfx4 file;
       dirIndex is the direction, from the cfx manual
       0 = solid (3-D patch),
       1 = high i, 2 = high j, 3 = high k
       4 = low i, 5 = low j, 6 = low k'''
    corrected_range = patch_range[:] # We won't change initial range
    
    for dirIndex in (dirIndex1, dirIndex2, dirIndex3):
        # high i
        if dirIndex == 1:   pass
        # high j
        elif dirIndex == 2: pass
        # high k
        elif dirIndex == 3: pass
        # low i
        elif dirIndex == 4:
            invert_indexes(corrected_range, 0, 1)
        # low j
        elif dirIndex == 5:
            invert_indexes(corrected_range, 2, 3)
        # low k
        elif dirIndex == 6:
            invert_indexes(corrected_range, 4, 5)
    
    return corrected_range
    

def main():
    logging.basicConfig(filename='cfx4ToPGSTK.log',
                        filemode='w',
                        format='%(asctime)s %(levelname)s: %(message)s',
                        level=logging.DEBUG)
        
    # Create console handler with a lower log level
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    # Create formatter and add it to the console handler
    formatter = logging.Formatter('%(asctime)s %(levelname)s: %(message)s')
    console_handler.setFormatter(formatter)
    # Add the handler to the logger
    logging.getLogger('').addHandler(console_handler)
    
    logger = logging.getLogger('cfx4ToPGSTK')
    
    if sys.argv[1]:
        cfxFile = open(sys.argv[1])
    else:
        logger.error("Please provide CFX4 file")
        sys.exit(3)
    
    logger.info("Reading PGSTK config")    
    
    pgs_config = PgstkConfig(prefix, working_dir, mesh_folder, output_folder)
    
    logger.info("Definition prefix: {prefix}".format(prefix=pgs_config.prefix))
    logger.info("Working directory: {wd}".format(wd = pgs_config.working_dir))
    logger.info("Mesh folder path: {mesh}".format(mesh = pgs_config.mesh_dir))
    logger.info("Output folder path: {out}".format(out = pgs_config.out_dir))
    
    logger.info("Reading blocks")
        
    if not cfxFile.readline().startswith("/*"):
        logger.error("First line should be a comment. Check file format")
        sys.exit(1)
    
    nblock, npatch, nglue, nelem, npoint = [int(_) for _ in cfxFile.readline().split()[:5]]
    logger.debug("NBLOCK = {}, NPATCH = {}, NGLUE = {}, NELEM = {}, NPOINT = {}".
                  format(nblock, npatch, nglue, nelem, npoint))
    
    if not cfxFile.readline().startswith("/*"):
        logger.error("Third line should be a comment. Check file format")
        sys.exit(1)
    
    hexBlock = {} # Structure for hexBlock is {blockNo_1:{blockName:"name1", NI:ni, NJ:nj, NK:nk},
                  #                            blockNo_2:{blockName:"name2", NI:ni, NJ:nj, NK:nk},
                  #                            ...
                  #                            blockNo_nblock:{blockName:"nameN", NI:ni, NJ:nj, NK:nk}}
    for blockNo in range(1, nblock+1):
        blockName, NI, NJ, NK = cfxFile.readline().split()
        hexBlock[blockNo] = dict(blockName=blockName, NI=int(NI), NJ=int(NJ), NK=int(NK))
    
    meshSize = 0
    for blockNo in list(hexBlock.keys()):
        NI, NJ, NK = (hexBlock[blockNo]["NI"],
                      hexBlock[blockNo]["NJ"],
                      hexBlock[blockNo]["NK"])
        hexBlock[blockNo]["blockCells"] = NI*NJ*NK
        hexBlock[blockNo]["blockVertices"] = (NI+1)*(NJ+1)*(NK+1)
        
        meshSize += hexBlock[blockNo]["blockCells"]
        
        if NI==1 or NJ==1 or NK==1:
            logger.error("Size of block {} is incorrect. Block number is {}. Check mesh".
                          format(hexBlock[blockNo]["blockName"], blockNo))
            sys.exit(2)

    logger.info("Mesh size is {:,} cells".format(meshSize))
    logger.info("Reasonable number of processes (NProc) for this task should be {:d}".format(int(meshSize/10000)))

    if not cfxFile.readline().startswith("/*"):
        logger.error("Line {} should be a comment. Check file format".format(cfxFile.tell()))
        sys.exit(1)
    
    logger.info("Reading patch definitions")
    
    cfxPatchType = []
    cfxPatchName = []
    cfxPatchNumber = []
    
    cfxPatchRange = []
    cfxPatchDirection = []
    cfxPatchBlkNo = []
    cfxPatchLabel = []
    
    eSlice = ExtendedSlice() # Instance for generating extended slices
    
    for patch in range(npatch):
        # Grab patch type, name and number
        patchType, patchName, patchNo = cfxFile.readline().split()
        cfxPatchType.append(patchType)
        cfxPatchName.append(patchName)
        cfxPatchNumber.append(int(patchNo))

        # Read next line and split it on slices
        patch_slice = get_multiple_slices(cfxFile.readline().split(), eSlice[:6,6:])
        
        # Grab patch range
        patchRange = [int(_) for _ in next(patch_slice)]
        cfxPatchRange.append(patchRange)

        # Grab direction, master block ID and patch label
        # Note: direc is the direction, from the cfx manual
        # 0 = solid (3-D patch),
        # 1 = high i, 2 = high j, 3 = high k
        # 4 = low i, 5 = low j, 6 = low k
        patchDirection, patchBlkNo, patchLabel = [int(_) for _ in next(patch_slice)]
        cfxPatchDirection.append(patchDirection)
        cfxPatchBlkNo.append(patchBlkNo)
        cfxPatchLabel.append(patchLabel)
        
    logger.debug("{} patch definitions readed".format(patch))
    
    if not cfxFile.readline().startswith("/*"):
        logger.error("Line {} should be a comment. Check file format".format(cfxFile.tell()))
        sys.exit(1)

    logger.info("Preparing folders for PGSTK definition")
    pgs_config.prepare_directories()
    pgs_config.set_blk_precision(nblock)
    pgs_config.set_default_bc(default_bc)

    logger.info("Reading block to block glueing information and creating cnn file")
    
    cnnFileFullPath = join(pgs_config.working_dir + os.sep, "{prefix}.cnn".format(prefix = pgs_config.prefix))
    cnnFile = open(cnnFileFullPath, 'w')
    
    glueMasterPatches = [] # We'll save glueing information in list structure just in case
    glueSlavePatches = []
    
    for glue in range(nglue):
        masterPatchNo, slavePatchNo, dirIndex1, dirIndex2, dirIndex3, joinNumber = [int(_) for _ in cfxFile.readline().split()]
        glueMasterPatches.append(masterPatchNo)
        glueSlavePatches.append(slavePatchNo)
        
        _index = lambda _: cfxPatchNumber.index(_) # Get index in cfxPatch array
        
        # Grab master patch range
        cnnMasterPatchRange = cells_to_vrt(cfxPatchRange[_index(masterPatchNo)],
                                           cfxPatchDirection[_index(masterPatchNo)])
        
        # Grab slave patch range
        cnnSlavePatchRange = cells_to_vrt(cfxPatchRange[_index(slavePatchNo)],
                                          cfxPatchDirection[_index(slavePatchNo)])
        # Correct patch diagonal
        cnnSlavePatchRange = correct_patch_range(cnnSlavePatchRange, dirIndex1, dirIndex2, dirIndex3)
        # Note: cnn file format is
        # 0                             ! 0 - matching grid; 1 - irregular connection
        # <block_name> <I1> <J1> <K1>   ! Start point
        #              <I2> <J2> <K2>   ! End point
        cnnFile.write("0\n")
        cnnFile.write("b{master_block:0{precision}d} {0} {2} {4} {1} {3} {5}\n".format(master_block = cfxPatchBlkNo[_index(masterPatchNo)],
                                                                                       precision=pgs_config.blk_precision,
                                                                                       *cnnMasterPatchRange))
        cnnFile.write("b{slave_block:0{precision}d} {0} {2} {4} {1} {3} {5}\n".format(slave_block = cfxPatchBlkNo[_index(slavePatchNo)],
                                                                                      precision=pgs_config.blk_precision,
                                                                                      *cnnSlavePatchRange))
    cnnFile.close()
    
    logger.debug("{} lines of glueing data readed".format(glue))
    
    logger.info("Reading block points and converting to PGSTK format")
    logger.info("Writing mesh files to appropriate folder and creating prefgrid.asg")
    
    prefgridFile = open(join(pgs_config.working_dir + os.sep, "prefgrid.asg"), 'w')
    prefgridFile.write("{prefix} {out} {mesh}\n".format(prefix = pgs_config.prefix,
                                                        out = basename(pgs_config.out_dir) + os.sep,
                                                        mesh = basename(pgs_config.mesh_dir) + os.sep))
    
    for blockNo in range(1, nblock+1):
        if not cfxFile.readline().startswith("/*"):
            logger.error("Line {} should be a comment. Check file format".format(cfxFile.tell()))
            sys.exit(1)

        meshFileName = 'b{:0{precision}d}.msh'.format(blockNo, precision=pgs_config.blk_precision)
        meshFileFullPath = join(pgs_config.mesh_dir + os.sep, meshFileName)
        meshFile = open(meshFileFullPath, 'w')
        
        prefgridFile.write("{grid_file} {block_name} {rank}\n".format(grid_file = meshFileName,
                                                                      block_name = splitext(basename(meshFileName))[0],
                                                                      rank = blockNo-1))
        meshFile.write("{} {} {}\n".format(hexBlock[blockNo]["NI"] + 1,
                                           hexBlock[blockNo]["NJ"] + 1,
                                           hexBlock[blockNo]["NK"] + 1))
        for point in range(hexBlock[blockNo]["blockVertices"]):
            meshFile.write("{} {} {}\n".format(*cfxFile.readline().split()))
        
        if not cfxFile.readline().startswith("/*"):
            logger.error("Line {} should be a comment. Check file format".format(cfxFile.tell()))
            sys.exit(1)
        
        meshFile.close()
    
    prefgridFile.close()

    cfxFile.close()
    
    logger.info("Searching for surface patches and writing them to container")
    
    volume_patches = "USER3D","POROUS","SOLID","SOLCON"
    connection_pathes = "BLKBDY",
    surface_patches = "WALL","SYMMET","INLET","OUTLET","USER2D"
    
    azn_list = [] # Structure for azn container: [(zone_name, zone_block, zone_range), ...]
    
    for patch in range(npatch):
        patch_type = cfxPatchType[patch]
        
        if patch_type.upper() in surface_patches:
            zone_name = cfxPatchName[patch]
            zone_block = cfxPatchBlkNo[patch]
            zone_range = cells_to_vrt(cfxPatchRange[patch],
                                      cfxPatchDirection[patch])
            
            azn_list.append((zone_name, zone_block, zone_range))
            
    azn_patches = defaultdict(list) # Using list as the default_factory
    
    # Generating structure for zones: {"Zone1_name":[(zone1_block1,zone1_range1),
    #                                                (zone1_block2,zone1_range2),...],
    #                                  ...
    #                                  "ZoneN_name":[(zoneN_block1,zoneN_range1),...]}    
    for zone_name, zone_block, zone_range in azn_list:
        azn_patches[zone_name].append((zone_block, zone_range))
    
    logger.info("Creating azn and bc files for current mesh")
    aznFileFullPath = join(pgs_config.working_dir + os.sep, "{prefix}.azn".format(prefix = pgs_config.prefix))
    bcFileFullPath = join(pgs_config.working_dir + os.sep, "{prefix}.bc".format(prefix = pgs_config.prefix))
    aznFile = open(aznFileFullPath, 'w')
    bcFile = open(bcFileFullPath, 'w')

    for zone in sorted(azn_patches.keys()):
        aznFile.write("%Zone% {zone_name} b{zone_block:0{precision}d} {0} {1} {2} {3} {4} {5}\n".format(zone_name = zone.lower(),
                                                                                                        zone_block = azn_patches[zone][0][0],
                                                                                                        precision=pgs_config.blk_precision,
                                                                                                        *azn_patches[zone][0][1]))
        bcFile.write("%Zone% {zone_name} {boundary_condition}\n".format(zone_name = zone.lower(),
                                                                        boundary_condition = pgs_config.default_bc))
        
        if len(azn_patches[zone]) > 1:
            for zone_block, zone_range in azn_patches[zone][1:]:
                aznFile.write("{space:{tab}} b{zone_block:0{precision}d} {0} {1} {2} {3} {4} {5}\n".format(space = " ",
                                                                                                           tab = len("%Zone% {}".format(zone)),
                                                                                                           zone_block = zone_block,
                                                                                                           precision=pgs_config.blk_precision,
                                                                                                           *zone_range))
    aznFile.close()
    bcFile.close()
    
    logger.info("Writing memory file")
    pgs_config.set_memory_vectors(RWK, IWK)
    memoryFileFullPath = join(pgs_config.working_dir + os.sep, "memory")
    memoryFile = open(memoryFileFullPath, 'w')
    memoryFile.write("{RWK} {IWK}\n".format(RWK = pgs_config.RWK, IWK =pgs_config.IWK))
    memoryFile.close()
    
    logger.info("Writing monitor file")
    monitorFileFullPath = join(pgs_config.working_dir + os.sep, "{prefix}.mnp".format(prefix = pgs_config.prefix))
    monitorFile = open(monitorFileFullPath, 'w')
    monitorFile.write("b{blk:0{precision}d} 0 0 0\n".format(blk=1, precision=pgs_config.blk_precision))
    monitorFile.close()
    
    logger.info("Mesh converted successfully")
    

if __name__ == '__main__':
    #Exit codes:
    #0 - success
    #1 - incorrect file format
    #2 - mesh is not compatible with PGSTK
    #3 - no input file specifyed
    main()