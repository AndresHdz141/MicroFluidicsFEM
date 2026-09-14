module setup



function setup_fevalues(ipu, ipp, ipg)
    qr = QuadratureRule{RefTriangle}(2)
    cvu = CellValues(qr, ipu, ipg)
    cvp = CellValues(qr, ipp, ipg; update_gradients = false, update_detJdV = false)
    qr_facet = FacetQuadratureRule{RefTriangle}(2)
    fvp = FacetValues(qr_facet, ipp, ipg; update_gradients = false)
    return cvu, cvp, fvp
end

end