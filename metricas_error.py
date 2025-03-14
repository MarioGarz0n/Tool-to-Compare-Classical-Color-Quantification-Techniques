import torch
import piq
from skimage.io import imread

import cv2 
import numpy as np

import sys 


# -------------------------------------------------------------
# Error cuadrático medio (MSE)
def mse (GT,P):
	"""calcula el error MSE (mean squared error).

	:param GT: primera imagen (original) 
	:param P: segunda imagen (modificada) 

	:returns:  float -- valor MSE.
	"""
	
	#GT,P = _initial_check(GT,P)


	return np.mean((GT.astype(np.float64)-P.astype(np.float64))**2)

# -------------------------------------------------------------
# Error absoluto medio
def mae (GT,P):
	"""calcula el error MAE (mean absolute error).

	:param GT: primera imagen (original) 
	:param P: segunda imagen (modificada)

	:returns:  float -- valor MAE.
	"""
	
	#GT,P = _initial_check(GT,P)


	return np.mean(np.abs(GT.astype(np.float64)-P.astype(np.float64)))
	
		
@torch.no_grad()
def main(figura1, figura2, N):
    # variable para indicar cómo se debe mostrar el resultado
    separar = 0 # 1 => separar los resultados en lineas diferentes
               # !=1 => mostrar todos los resultados en la misma linea
               
               
    # Se leen las dos imágenes RGB
    x = torch.tensor(imread(figura1)).permute(2, 0, 1)[None, ...] 
    y = torch.tensor(imread(figura2)).permute(2, 0, 1)[None, ...] 
    

    if torch.cuda.is_available():
        # Move to GPU to make computaions faster
        x = x.cuda()
        y = y.cuda()

    #print(x.size())  IMPRIME: torch.Size([1, 3, 512, 512])


    # se calculan varios índices de error
    # -- To compute BRISQUE score as a measure, use lower case function from the library
    #brisque_index: torch.Tensor = piq.brisque(x, data_range=1., reduction='none')
    brisque_index: torch.Tensor = piq.brisque(x, data_range=255, reduction='none')    
 
    # -- DSS (discrete cosine transform Subbands Similarity>)
#    dss_index: torch.Tensor = piq.dss(x, y, data_range=1., reduction='none')
    dss_index: torch.Tensor = piq.dss(x, y, data_range=255, reduction='none')    

    # -- FSIM 
    fsim_index: torch.Tensor = piq.fsim(x, y, data_range=255, reduction='none')

    # -- To compute GMSD as a measure, use lower case function from the library
    # This is port of MATLAB version from the authors of original paper.
    # In any case it should me minimized. Usually values of GMSD lie in [0, 0.35] interval.
    gmsd_index: torch.Tensor = piq.gmsd(x, y, data_range=255, reduction='none')

    # -- HaarPSI 
    # This is port of MATLAB version from the authors of original paper.
    haarpsi_index: torch.Tensor = piq.haarpsi(x, y, data_range=255, reduction='none')

    # -- MDSI 
    mdsi_index: torch.Tensor = piq.mdsi(x, y, data_range=255, reduction='none')
    
    # -- To compute Multi-Scale GMSD as a measure, use lower case function from the library
    # It can be used both as a measure and as a loss function. In any case it should Be minimized.
    # By default scale weights are initialized with values from the paper.
    # You can change them by passing a list of 4 variables to scale_weights argument during initialization
    # Note that input tensors should contain images with height and width equal 2 ** number_of_scales + 1 at least.
    ms_gmsd_index: torch.Tensor = piq.multi_scale_gmsd(
        x, y, data_range=255, chromatic=True, reduction='none')
    
    # -- PSNR 
    psnr_index = piq.psnr(x, y, data_range=255, reduction='none')
    
    # -- SSIM index 
    ssim_index = piq.ssim(x, y,  data_range=255, kernel_size = 11) #, k1 = 0.01, k2 = 0.02)    

    # -- MS-SSIM index 
    ms_ssim_index: torch.Tensor = piq.multi_scale_ssim(x, y, data_range=255)

    # -- W-SSIM index 
    wissim_index = piq.information_weighted_ssim(x, y, data_range=255)

    # -- Error MSE, cuyo calculo se ha añadido a las funciones de la librería
    mse_index = piq.mse(x, y, data_range=255, reduction='none')

    # -- TV 
    tv_index: torch.Tensor = piq.total_variation(x)

    # -- VIF 
    vif_index: torch.Tensor = piq.vif_p(x, y, data_range=255)

    # -- VSI 
    vsi_index: torch.Tensor = piq.vsi(x, y, data_range=255)

    # -- SR-SIM score 
    srsim_index: torch.Tensor = piq.srsim(x, y, data_range=255)


    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    #
    # To compute Style score as a loss function, use corresponding PyTorch module:
    # By default VGG16 model is used, but any feature extractor model is supported.
    # Don't forget to adjust layers names accordingly. Features from different layers can be weighted differently.
    # Use weights parameter. See other options in class docstring.
    #OJO    QUITO ESTO, QUE ME INTENTA DESCARGAR EL MODELO VGG16    
#    style_loss = piq.StyleLoss(feature_extractor="vgg16", layers=("relu3_3",))(x, y)
#    print(f"Style: {style_loss.item():0.4f}")


    # To compute Content score as a loss function, use corresponding PyTorch module
    # By default VGG16 model is used, but any feature extractor model is supported.
    # Don't forget to adjust layers names accordingly. Features from different layers can be weighted differently.
    # Use weights parameter. See other options in class docstring.
    #OJO    QUITO ESTO, QUE ME INTENTA DESCARGAR EL MODELO VGG16
    #content_loss = piq.ContentLoss(
    #    feature_extractor="vgg16", layers=("relu3_3",), reduction='none')(x, y)
    #print(f"ContentLoss: {content_loss.item():0.4f}")


    # To compute DISTS as a loss function, use corresponding PyTorch module
    # By default input images are normalized with ImageNet statistics before forwarding through VGG16 model.
    # If there is no need to normalize the data, use mean=[0.0, 0.0, 0.0] and std=[1.0, 1.0, 1.0].
    #dists_loss = piq.DISTS(reduction='none')(x, y)
    #print(f"DISTS: {dists_loss.item():0.4f}")


    # To compute LPIPS as a loss function, use corresponding PyTorch module
    #lpips_loss: torch.Tensor = piq.LPIPS(reduction='none')(x, y)
    #print(f"LPIPS: {lpips_loss.item():0.4f}")
  
  
    # To compute PieAPP as a loss function, use corresponding PyTorch module:
    #pieapp_loss: torch.Tensor = piq.PieAPP(reduction='none', stride=32)(x, y)
    #print(f"PieAPP loss: {pieapp_loss.item():0.4f}")



    FIG1 = cv2.imread(figura1,cv2.IMREAD_COLOR)  
    FIG2 = cv2.imread(figura2,cv2.IMREAD_COLOR)   

    resul_mse= mse (FIG1, FIG2)
    
    mse_res = resul_mse * 3

    resul_mae = 3*mae(FIG1, FIG2)


    # se muestran los índices
    if separar == 1:
       print(figura1 + " "+ figura2) # nombre de ambas imágenes
       
       print ('MSE:', f"{mse_res:0.5f}") 
       print(f"PSNR: {psnr_index.item():0.5f}")        
       print ('MAE:', f"{resul_mae:05f}")    
                
       #print(f"BRISQUE: {brisque_index.item():0.5f}")
       print(f"DSS: {dss_index.item():0.5f}")
       print(f"FSIM: {fsim_index.item():0.5f}")
       print(f"GMSD: {gmsd_index.item():0.5f}")
       print(f"HaarPSI: {haarpsi_index.item():0.5f}")
       print(f"MDSI: {mdsi_index.item():0.5f}")
       print(f"MS-GMSDc: {ms_gmsd_index.item():0.5f}")    
       print(f"SSIM: {ssim_index.item():0.5f}")    
       print(f"MS-SSIM: {ms_ssim_index.item():0.5f}")    
       print(f"iw-SSIM: {wissim_index.item():0.5f}")    
       #print(f"TV index: {tv_index.item():0.5f}")
       print(f"VIFp: {vif_index.item():0.5f}")    
       print(f"VSI: {vsi_index.item():0.5f}")    
       print(f"SR-SIM: {srsim_index.item():0.5f}")    
       
    else:
       # pongo la palabra 'imagen' para poder alinear los errores con sus etiquetas
       # al copiar los resultados a una hoja de cálculo
       #print("imagen ", end=" ")
              
       # muestro las etiquetas
       #print ('MSE', end=" ") 
       #print(f"PSNR ", end=" ")          
       #print ('MAE', end=" ")  
              
       #print("BRISQUE ", end=" ")
       #print("DSS ",end=" ")
       #print("FSIM ", end=" ")
       #print(f"GMSD ", end=" ")
       #print(f"HaarPSI ", end=" ")
       #print(f"MDSI ", end=" ")
       #print(f"MS-GMSDc ", end=" ")    
       #print(f"MSE ", end=" ")    
       #print(f"SSIM ", end=" ")    
       #print(f"MS-SSIM ", end=" ")    
       #print(f"iw-SSIM ", end=" ")    
       #print(f"TV index ", end=" ")
       #print(f"VIFp ", end=" ")    
       #print(f"VSI ", end=" ")    
       #print(f"SR-SIM ")    

       print(N, end=" ") #numero de colores de la paleta
       print(figura1, end=" ")  #nombre de la primera de las imágenes       
       
       #muestro los valores
       print ( f"{mse_res:0.5f}", end=" ") 
       print(f"{psnr_index.item():0.5f}", end=" ")        
       print ( f"{resul_mae:05f}", end=" ")  
       
       #print(f"{brisque_index.item():0.5f}", end=" ")
       print(f"{dss_index.item():0.5f}", end=" ")
       print(f"{fsim_index.item():0.5f}", end=" ")
       print(f"{gmsd_index.item():0.5f}", end=" ")
       print(f"{haarpsi_index.item():0.5f}", end=" ")
       print(f"{mdsi_index.item():0.5f}", end=" ")
       print(f"{ms_gmsd_index.item():0.5f}", end=" ")  
       print(f"{ssim_index.item():0.5f}", end=" ")    
       print(f"{ms_ssim_index.item():0.5f}", end=" ")    
       print(f"{wissim_index.item():0.5f}", end=" ")    
       #print(f"{tv_index.item():0.5f}", end=" ")
       print(f"{vif_index.item():0.5f}", end=" ")    
       print(f"{vsi_index.item():0.5f}", end=" ")    
       print(f"{srsim_index.item():0.5f}")    



if __name__ == '__main__':
    #main()
    if len(sys.argv) == 4:
      figura1 = sys.argv[1]
      figura2 = sys.argv[2]
      N = sys.argv[3]

      #print("procesando las imagenes: |"+ figura1 + "| y |"+ figura2 + "|")
      #print(figura1 + " "+ figura2) #NOMBRES DE LOS DOS FICHEROS

      
      main(figura1, figura2, N)
    else:
      print("Debes indicar la imagen original y la cuantizada")
