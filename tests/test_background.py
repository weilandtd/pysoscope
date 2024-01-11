import scanpy as sc

from pysoscope.analysis.background import subtract_background, find_background_pixels, find_background_signature


def test_background():
    adata = sc.read_h5ad("test_anndata.h5ad")

    # Test if number on picels is plaubsible
    background_pixels = find_background_pixels(adata, [124.007,], 25)
    print(f"{background_pixels.sum()} out of {len(background_pixels)} pixels are background pixels")

    # Test if background signature is plausible
    background_signature = find_background_signature(adata, background_pixels, 25)
    print(background_signature)

    # This is the only requiered function for background substraction
    adata_mod = subtract_background(adata, [124.007,], 75, 25)

    