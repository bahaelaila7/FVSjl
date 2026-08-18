# LPOPDY sub-chunk 1: MPBDRV phloem model — bit-exact vs FVSie_lpmpb_ph (fort.779 dump).
# DDS5=(2*DBH*DG+DG^2)/2; BAI5=DDS5*0.7853982; if BAI5>1e-4:
#   XPT=exp(-3.17152 + .12591*log(BAI5) + .50932*log(DBH) - .0077*HT)   (glibc expf/logf)
h2f(h)=reinterpret(Float32,parse(UInt32,h;base=16))
f2h(x)=uppercase(string(reinterpret(UInt32,Float32(x));base=16,pad=8))
lf(x)=ccall(:logf,Float32,(Float32,),x); ef(x)=ccall(:expf,Float32,(Float32,),x)
function mpb_phloem(dbh::Float32,dg::Float32,ht::Float32)::Float32
    dds5=(2f0*dbh*dg + dg*dg)/2f0; bai5=dds5*0.7853982f0
    bai5 <= 0.0001f0 && return 0f0
    return ef(-3.17152f0 + 0.12591f0*lf(bai5) + 0.50932f0*lf(dbh) - 0.0077f0*ht)
end
function main()
  ok=0; n=0
  for ln in eachline("/workspace/.iework/lpmpb/run/fort.779")
    f=split(strip(ln)); length(f)<5 && continue; n+=1
    dbh,dg,ht,_=h2f.(f[2:5]); ok += (f2h(mpb_phloem(dbh,dg,ht))==f[5])
  end
  println("LPOPDY phloem: $ok/$n bit-exact")
end
main()
