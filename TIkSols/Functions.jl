using PyPlot
function reduced_diff(vdm,tik)
  sum(abs(vdm-tik))/length(vdm)
end


function disp(vdm,tik)
    figure()
    subplot(131)
    imshow(vdm,interpolation="none",cmap="Greys")
    title("VDM")

    subplot(132)
    imshow(tik,interpolation="none",cmap="Greys")
    title("TIK")

    subplot(133)
    imshow(abs(tik-vdm),interpolation="none",cmap="Greys")
    title("Difference")

    draw()
  end

function evaldisp(res)
  mu=res[:,1]
  chi2=res[:,2]
  diff=res[:,3]
  figure()
  subplot(2,1,1)
  loglog(mu,chi2,".")
  xlabel("Mu")
  ylabel("Chi2")
  subplot(2,1,2)
  semilogx(mu,diff,".")
  xlabel("Mu")
  ylabel("Diff")
  tight_layout()
  draw()
end
