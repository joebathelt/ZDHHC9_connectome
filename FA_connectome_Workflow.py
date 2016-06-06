#=====================================================================================
# Additional nodes for the workflow

from nipype.interfaces.base import BaseInterface, BaseInterfaceInputSpec, CommandLineInputSpec, CommandLine, traits, File, TraitedSpec
from nipype.interfaces.matlab import MatlabCommand

#==================================================================================================
# Denoising with non-local means
# This function is based on the example in the Dipy preprocessing tutorial:
# http://nipy.org/dipy/examples_built/denoise_nlmeans.html#example-denoise-nlmeans

class DipyDenoiseInputSpec(BaseInterfaceInputSpec):
	in_file = File(exists=True, desc='diffusion weighted volume for denoising', mandatory=True)

class DipyDenoiseOutputSpec(TraitedSpec):
	out_file = File(exists=True, desc="denoised diffusion-weighted volume")

class DipyDenoise(BaseInterface):
	input_spec = DipyDenoiseInputSpec
	output_spec = DipyDenoiseOutputSpec

	def _run_interface(self, runtime):
		import nibabel as nib
		import numpy as np
		import matplotlib.pyplot as plt
		from dipy.denoise.nlmeans import nlmeans
		from nipype.utils.filemanip import split_filename

		fname = self.inputs.in_file
		img = nib.load(fname)
		data = img.get_data()
		affine = img.get_affine()
		mask = data[..., 0] > 80
		a = data.shape 

		denoised_data = np.ndarray(shape=data.shape)
		for image in range(0,a[3]):
		    print(str(image + 1) + '/' + str(a[3] + 1))
		    dat = data[...,image]
		    sigma = np.std(dat[~mask]) # Calculating the standard deviation of the noise 
		    den = nlmeans(dat, sigma=sigma, mask=mask)
		    denoised_data[:,:,:,image] = den

		_, base, _ = split_filename(fname)
		nib.save(nib.Nifti1Image(denoised_data, affine), base + '_denoised.nii')

		return runtime

	def _list_outputs(self):
		from nipype.utils.filemanip import split_filename
		import os 
		outputs = self._outputs().get()
		fname = self.inputs.in_file
		_, base, _ = split_filename(fname)
		outputs["out_file"] = os.path.abspath(base + '_denoised.nii')
		return outputs

#==================================================================================================
# Moving tracts to a different space
class trk_CoregInputSpec(CommandLineInputSpec):
	in_file = File(exists=True, desc='whole-brain tractography in .trk format', 
		mandatory=True, position = 0, argstr="%s")
	output_file = File("coreg_tracks.trk", desc="whole-brain tractography in coregistered space", 
		position=1, argstr="%s", usedefault=True)
	FA_file = File(exists=True, desc='FA file in the same space as the .trk file', 
		mandatory=True, position = 2, argstr="-src %s")
	reference = File(exists=True, desc='Image that the .trk file will be registered to', 
		mandatory=True, position = 3, argstr="-ref %s")
	transfomation_matrix = File(exists=True, desc='FSL matrix with transform form original to new space', 
		mandatory=True, position = 4, argstr="-reg %s")

class trk_CoregOutputSpec(TraitedSpec):
	transformed_track_file = File(exists=True, desc="whole-brain tractography in new space")

class trk_Coreg(CommandLine):
	input_spec = trk_CoregInputSpec
	output_spec = trk_CoregOutputSpec

	_cmd = "track_transform"

	def _list_outputs(self):#
		import os 
		outputs = self.output_spec().get()
		outputs['transformed_track_file'] = os.path.abspath(self.inputs.output_file)
		return outputs

#==================================================================================================
# Extract b0
class Extractb0InputSpec(BaseInterfaceInputSpec):
	in_file = File(exists=True, desc='diffusion-weighted image (4D)', mandatory=True)

class Extractb0OutputSpec(TraitedSpec):
	out_file = File(exists=True, desc="First volume of the dwi file")

class Extractb0(BaseInterface):
	input_spec = Extractb0InputSpec
	output_spec = Extractb0OutputSpec

	def _run_interface(self, runtime):
		import nibabel as nib
		img = nib.load(self.inputs.in_file)
		data = img.get_data()
		affine = img.get_affine()

		from nipype.utils.filemanip import split_filename
		import os 
		outputs = self._outputs().get()
		fname = self.inputs.in_file
		_, base, _ = split_filename(fname)
		nib.save(nib.Nifti1Image(data[...,0],affine),os.path.abspath(base + '_b0.nii.gz'))
		return runtime

	def _list_outputs(self):
		from nipype.utils.filemanip import split_filename
		import os 
		outputs = self._outputs().get()
		fname = self.inputs.in_file
		_, base, _ = split_filename(fname)
		outputs["out_file"] = os.path.abspath(base + '_b0.nii.gz')
		return outputs
			
#==================================================================================================
# FA connectome construction

class FAconnectomeInputSpec(BaseInterfaceInputSpec):
	trackfile = File(exists=True, desc='whole-brain tractography in .trk format', mandatory=True)
	ROI_file = File(exists=True, desc='image containing the ROIs', mandatory=True)
	FA_file = File(exists=True, desc='fractional anisotropy map in the same soace as the track file', mandatory=True)
	output_file = File("FA_matrix.txt", desc="Adjacency matrix of ROIs with FA as conenction weight", usedefault=True)

class FAconnectomeOutputSpec(TraitedSpec):
	out_file = File(exists=True, desc="connectivity matrix of FA between each pair of ROIs")

class FAconnectome(BaseInterface):
	input_spec = FAconnectomeInputSpec
	output_spec = FAconnectomeOutputSpec

	def _run_interface(self, runtime):
		# Loading the ROI file
	    import nibabel as nib
	    import numpy as np
	    from dipy.tracking import utils 

	    img = nib.load(self.inputs.ROI_file)
	    data = img.get_data()
	    affine = img.get_affine()

	    # Getting ROI volumes if they haven't been generated
		if not os.path.isfile('/imaging/jb07/CALM/DWI/FA_connectome/Atlas_volumes.csv'):
			import nibabel as nib
			import numpy as np
			import os 
			import pandas as pd
			import subprocess

			atlas_file = ROI_file
			img = nib.load(atlas_file)
			data = img.get_data()
			affine = img.get_affine()
			volumes = pd.DataFrame()

			atlas_labels = np.unique(data)

			for atlas_label in atlas_labels:
				data = nib.load(atlas_file).get_data()
				data[data != atlas_label] = 0 
				data[data == atlas_label] = 1
				nib.save(nib.Nifti1Image(data, affine), 'temp.nii.gz')
				volumes.set_value(atlas_label, 'volume', subprocess.check_output(os.environ['FSLDIR'] + '/bin/fslstats temp.nii.gz -V', shell=True).split(' ')[0])

			os.remove('temp.nii.gz')
			volumes.to_csv('/imaging/jb07/CALM/DWI/FA_connectome/Atlas_volumes.csv')

		ROI_volumes = pd.read_csv('/home/jb07/CALM/DWI/FA_connectome/Atlas_volumes.csv')

		# Getting the FA file
		img = nib.load(FA_file)
		FA_data = img.get_data()
		FA_affine = img.get_affine()

		# Loading the streamlines
		from nibabel import trackvis
		streams, hdr = trackvis.read(trackfile,points_space='rasmm')
		streamlines = [s[0] for s in streams]
		streamlines_affine = trackvis.aff_from_hdr(hdr,atleast_v2=True)

		# Checking for negative values
		from dipy.tracking._utils import _mapping_to_voxel, _to_voxel_coordinates
		endpoints = [sl[0::len(sl)-1] for sl in streamlines]
		lin_T, offset = _mapping_to_voxel(affine, (1.,1.,1.))
		inds = np.dot(endpoints, lin_T)
		inds += offset
		negative_values = np.where(inds <0)[0]
		for negative_value in sorted(negative_values, reverse=True):
			del streamlines[negative_value]

		# Constructing the streamlines matrix
		matrix,mapping = utils.connectivity_matrix(streamlines=streamlines,label_volume=data,affine=streamlines_affine,symmetric=True,return_mapping=True,mapping_as_streamlines=True)
		matrix[matrix < 10] = 0

		# Constructing the FA matrix
		dimensions = matrix.shape
		FA_matrix = np.empty(shape=dimensions)
		density_matrix = np.empty(shape=dimensions)
		density_corrected_matrix = np.empty(shape=dimensions)

		for i in range(0,dimensions[0]):
		    for j in range(0,dimensions[1]):
		        if matrix[i,j]:
		            dm = utils.density_map(mapping[i,j], FA_data.shape, affine=streamlines_affine)
		            FA_matrix[i,j] = np.mean(FA_data[dm>0])
		            if np.sum(dm > 0) > 0:
		            	density_matrix[i,j] = np.sum(dm[dm > 0])
		            	density_corrected_matrix[i,j] = np.sum(dm[dm > 0])/np.sum([ROI_volumes.iloc[i].values.astype('int'), ROI_volumes.iloc[j].values.astype('int')])
		            else: 
		            	density_matrix[i,j] = 0
		            	density_corrected_matrix[i,j] = 0
		        else:
		            FA_matrix[i,j] = 0
		            density_matrix[i,j] = 0
		            density_corrected_matrix[i,j] = 0

		FA_matrix[np.tril_indices(n=len(FA_matrix))] = 0
		FA_matrix = FA_matrix.T + FA_matrix - np.diagonal(FA_matrix)

		density_matrix[np.tril_indices(n=len(density_matrix))] = 0
		density_matrix = density_matrix.T + density_matrix - np.diagonal(density_matrix)

		density_corrected_matrix[np.tril_indices(n=len(density_corrected_matrix))] = 0
		density_corrected_matrix = density_corrected_matrix.T + density_corrected_matrix - np.diagonal(density_corrected_matrix)

		from nipype.utils.filemanip import split_filename
		_, base, _ = split_filename(self.inputs.trackfile)
		np.savetxt(base + '_FA_matrix.txt',FA_matrix,delimiter='\t')
		np.savetxt(base + '_density_matrix.txt',density_matrix,delimiter='\t')
		np.savetxt(base + '_volume_corrected_density_matrix.txt',density_corrected_matrix,delimiter='\t')
	    return runtime

	def _list_outputs(self):
		from nipype.utils.filemanip import split_filename
		import os 
		outputs = self._outputs().get()
		fname = self.inputs.trackfile
		_, base, _ = split_filename(fname)
		outputs["out_file"] = os.path.abspath(base + '_FA_matrix.txt')
		return outputs


#==================================================================================================
# Convert an adjacency matrix in txt format to NetworkX pck format

class TXT2PCKInputSpec(BaseInterfaceInputSpec):
	in_file = File(exists=True, desc='adjacency matrix in txt format', mandatory=True)

class TXT2PCKOutputSpec(TraitedSpec):
	out_file = File(exists=True, desc="NetworkX file in pck format")

class TXT2PCK(BaseInterface):
	input_spec = TXT2PCKInputSpec
	output_spec = TXT2PCKOutputSpec

	def _run_interface(self, runtime):
		# Reading the matrix file
		import numpy as np
		import networkx as nx

		adjacency_matrix = np.loadtxt(self.inputs.in_file)
		G = nx.from_numpy_matrix(adjacency_matrix)

		from nipype.utils.filemanip import split_filename
		_, base, _ = split_filename(self.inputs.in_file)
		nx.write_gpickle(G,path=base + '.pck')
		return runtime

	def _list_outputs(self):
		from nipype.utils.filemanip import split_filename
		import os 
		outputs = self._outputs().get()
		fname = self.inputs.in_file
		_, base, _ = split_filename(fname)
		outputs["out_file"] = os.path.abspath(base + '.pck')
		return outputs

#===============================================================================
# FA connectome workflow 

def FA_connectome(subject_list,base_directory,out_directory):

	#==============================================================
	# Loading required packages
	import nipype.interfaces.io as nio
	import nipype.pipeline.engine as pe
	import nipype.interfaces.utility as util
	import nipype.interfaces.fsl as fsl
	import nipype.interfaces.dipy as dipy
	import nipype.interfaces.mrtrix as mrt
	from own_nipype import DipyDenoise as denoise
	from own_nipype import trk_Coreg as trkcoreg
	from own_nipype import TXT2PCK as txt2pck
	from own_nipype import FAconnectome as connectome
	from own_nipype import Extractb0 as extract_b0
	import nipype.interfaces.cmtk as cmtk
	import nipype.interfaces.diffusion_toolkit as dtk
	import nipype.algorithms.misc as misc

	from nipype import SelectFiles
	import os
	registration_reference = os.environ['FSLDIR'] + '/data/standard/MNI152_T1_1mm_brain.nii.gz'
	nodes = list()

	#====================================
	# Defining the nodes for the workflow

	# Utility nodes
	gunzip = pe.Node(interface=misc.Gunzip(), name='gunzip')
	gunzip2 = pe.Node(interface=misc.Gunzip(), name='gunzip2')
	fsl2mrtrix = pe.Node(interface=mrt.FSL2MRTrix(invert_x=True),name='fsl2mrtrix')

	# Getting the subject ID
	infosource  = pe.Node(interface=util.IdentityInterface(fields=['subject_id']),name='infosource')
	infosource.iterables = ('subject_id', subject_list)

	# Getting the relevant diffusion-weighted data
	templates = dict(dwi='{subject_id}/dwi/{subject_id}_dwi.nii.gz',
		bvec='{subject_id}/dwi/{subject_id}_dwi.bvec',
		bval='{subject_id}/dwi/{subject_id}_dwi.bval')

	selectfiles = pe.Node(SelectFiles(templates),
	                   name='selectfiles')
	selectfiles.inputs.base_directory = os.path.abspath(base_directory)

	# Denoising
	denoise = pe.Node(interface=denoise(), name='denoise')

	# Eddy-current and motion correction
	eddycorrect = pe.Node(interface=fsl.epi.EddyCorrect(), name='eddycorrect')
	eddycorrect.inputs.ref_num = 0

	# Upsampling
	resample = pe.Node(interface=dipy.Resample(interp=3,vox_size=(1.,1.,1.)), name='resample')

	# Extract b0 image
	extract_b0 = pe.Node(interface=extract_b0(),name='extract_b0')

	# Fitting the diffusion tensor model
	dwi2tensor = pe.Node(interface=mrt.DWI2Tensor(), name='dwi2tensor')
	tensor2vector = pe.Node(interface=mrt.Tensor2Vector(), name='tensor2vector')
	tensor2adc = pe.Node(interface=mrt.Tensor2ApparentDiffusion(), name='tensor2adc')
	tensor2fa = pe.Node(interface=mrt.Tensor2FractionalAnisotropy(), name='tensor2fa')

	# Create a brain mask
	bet = pe.Node(interface=fsl.BET(frac=0.3,robust=False,mask=True),name='bet')

	# Eroding the brain mask
	erode_mask_firstpass = pe.Node(interface=mrt.Erode(), name='erode_mask_firstpass')
	erode_mask_secondpass = pe.Node(interface=mrt.Erode(), name='erode_mask_secondpass')
	MRmultiply = pe.Node(interface=mrt.MRMultiply(), name='MRmultiply')
	MRmult_merge = pe.Node(interface=util.Merge(2), name='MRmultiply_merge')
	threshold_FA = pe.Node(interface=mrt.Threshold(absolute_threshold_value = 0.7), name='threshold_FA')

	# White matter mask
	gen_WM_mask = pe.Node(interface=mrt.GenerateWhiteMatterMask(), name='gen_WM_mask')
	threshold_wmmask = pe.Node(interface=mrt.Threshold(absolute_threshold_value = 0.4), name='threshold_wmmask')

	# CSD probabilistic tractography 
	estimateresponse = pe.Node(interface=mrt.EstimateResponseForSH(maximum_harmonic_order = 8), name='estimateresponse')
	csdeconv = pe.Node(interface=mrt.ConstrainedSphericalDeconvolution(maximum_harmonic_order = 8), name='csdeconv')

	# Tracking 
	probCSDstreamtrack = pe.Node(interface=mrt.ProbabilisticSphericallyDeconvolutedStreamlineTrack(), name='probCSDstreamtrack')
	probCSDstreamtrack.inputs.inputmodel = 'SD_PROB'
	probCSDstreamtrack.inputs.desired_number_of_tracks = 150000
	tck2trk = pe.Node(interface=mrt.MRTrix2TrackVis(), name='tck2trk')

	# smoothing the tracts 
	smooth = pe.Node(interface=dtk.SplineFilter(step_length=0.5), name='smooth')

	# Co-registration with MNI space
	mrconvert = pe.Node(mrt.MRConvert(extension='nii'), name='mrconvert')
	flt = pe.Node(interface=fsl.FLIRT(reference=registration_reference, dof=12, cost_func='corratio'), name='flt')

	# Moving tracts to common space
	trkcoreg = pe.Node(interface=trkcoreg(reference=registration_reference),name='trkcoreg')

	# calcuating the connectome matrix 
	calc_matrix = pe.Node(interface=connectome(ROI_file='/home/jb07/Desktop/aal.nii.gz'),name='calc_matrix')

	# Converting the adjacency matrix from txt to pck format
	txt2pck = pe.Node(interface=txt2pck(), name='txt2pck')

	# Calculate graph theory measures with NetworkX and CMTK
	nxmetrics = pe.Node(interface=cmtk.NetworkXMetrics(treat_as_weighted_graph = True), name='nxmetrics')

	#====================================
	# Setting up the workflow
	fa_connectome = pe.Workflow(name='FA_connectome')

	# Reading in files
	fa_connectome.connect(infosource, 'subject_id', selectfiles, 'subject_id')

	# Denoising
	fa_connectome.connect(selectfiles, 'dwi', denoise, 'in_file')

	# Eddy current and motion correction
	fa_connectome.connect(denoise, 'out_file',eddycorrect, 'in_file')
	fa_connectome.connect(eddycorrect, 'eddy_corrected', resample, 'in_file')
	fa_connectome.connect(resample, 'out_file', extract_b0, 'in_file')
	fa_connectome.connect(resample, 'out_file', gunzip,'in_file')

	# Brain extraction
	fa_connectome.connect(extract_b0, 'out_file', bet, 'in_file')

	# Creating tensor maps
	fa_connectome.connect(selectfiles,'bval',fsl2mrtrix,'bval_file')
	fa_connectome.connect(selectfiles,'bvec',fsl2mrtrix,'bvec_file')
	fa_connectome.connect(gunzip,'out_file',dwi2tensor,'in_file')
	fa_connectome.connect(fsl2mrtrix,'encoding_file',dwi2tensor,'encoding_file')
	fa_connectome.connect(dwi2tensor,'tensor',tensor2vector,'in_file')
	fa_connectome.connect(dwi2tensor,'tensor',tensor2adc,'in_file')
	fa_connectome.connect(dwi2tensor,'tensor',tensor2fa,'in_file')
	fa_connectome.connect(tensor2fa,'FA', MRmult_merge, 'in1')

	# Thresholding to create a mask of single fibre voxels
	fa_connectome.connect(gunzip2, 'out_file', erode_mask_firstpass, 'in_file')
	fa_connectome.connect(erode_mask_firstpass, 'out_file', erode_mask_secondpass, 'in_file')
	fa_connectome.connect(erode_mask_secondpass,'out_file', MRmult_merge, 'in2')
	fa_connectome.connect(MRmult_merge, 'out', MRmultiply,  'in_files')
	fa_connectome.connect(MRmultiply, 'out_file', threshold_FA, 'in_file')

	# Create seed mask
	fa_connectome.connect(gunzip, 'out_file', gen_WM_mask, 'in_file')
	fa_connectome.connect(bet, 'mask_file', gunzip2, 'in_file')
	fa_connectome.connect(gunzip2, 'out_file', gen_WM_mask, 'binary_mask')
	fa_connectome.connect(fsl2mrtrix, 'encoding_file', gen_WM_mask, 'encoding_file')
	fa_connectome.connect(gen_WM_mask, 'WMprobabilitymap', threshold_wmmask, 'in_file')

	# Estimate response
	fa_connectome.connect(gunzip, 'out_file', estimateresponse, 'in_file')
	fa_connectome.connect(fsl2mrtrix, 'encoding_file', estimateresponse, 'encoding_file')
	fa_connectome.connect(threshold_FA, 'out_file', estimateresponse, 'mask_image')

	# CSD calculation
	fa_connectome.connect(gunzip, 'out_file', csdeconv, 'in_file')
	fa_connectome.connect(gen_WM_mask, 'WMprobabilitymap', csdeconv, 'mask_image')
	fa_connectome.connect(estimateresponse, 'response', csdeconv, 'response_file')
	fa_connectome.connect(fsl2mrtrix, 'encoding_file', csdeconv, 'encoding_file')

	# Running the tractography
	fa_connectome.connect(threshold_wmmask, "out_file", probCSDstreamtrack, "seed_file")
	fa_connectome.connect(csdeconv, "spherical_harmonics_image", probCSDstreamtrack, "in_file")
	fa_connectome.connect(gunzip, "out_file", tck2trk, "image_file")
	fa_connectome.connect(probCSDstreamtrack, "tracked", tck2trk, "in_file")

	# Smoothing the trackfile
	fa_connectome.connect(tck2trk, 'out_file',smooth,'track_file')

	# Co-registering FA with FMRIB58_FA_1mm standard space 
	fa_connectome.connect(MRmultiply,'out_file',mrconvert,'in_file')
	fa_connectome.connect(mrconvert,'converted',flt,'in_file')
	fa_connectome.connect(smooth,'smoothed_track_file',trkcoreg,'in_file')
	fa_connectome.connect(mrconvert,'converted',trkcoreg,'FA_file')
	fa_connectome.connect(flt,'out_matrix_file',trkcoreg,'transfomation_matrix')

	# Calculating the FA connectome
	fa_connectome.connect(trkcoreg,'transformed_track_file',calc_matrix,'trackfile')
	fa_connectome.connect(flt,'out_file',calc_matrix,'FA_file')

	# Calculating graph measures 
	fa_connectome.connect(calc_matrix,'out_file',txt2pck,'in_file')
	fa_connectome.connect(txt2pck,'out_file',nxmetrics,'in_file')

	#====================================
	# Running the workflow
	fa_connectome.base_dir = os.path.abspath(out_directory)
	fa_connectome.write_graph()
	fa_connectome.run('PBSGraph')

def get_folder_names(folder):
	import os,re 
	folder_names = list()

	for subfolder in os.listdir(folder):
		if re.search('CBU',subfolder):
			folder_names.append(subfolder)

	return folder_names