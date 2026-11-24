﻿﻿﻿﻿# Number audit: every printed statistic must agree with the archived tables at the
# precision printed. Stops the build instead of shipping a stale number.
Push-Location 'D:\ar_pilot'
python 'D:/ar_pilot/build/audit_numbers.py' --list | Tee-Object -FilePath 'D:/ar_pilot/build/audit_last.txt'
if ($LASTEXITCODE -ne 0) { Pop-Location; throw "audit_numbers.py reported a mismatch - build stopped" }
Pop-Location
# Rebuild both PDFs from source.
# Needs: pandoc (winget install JohnMacFarlane.Pandoc) and MiKTeX (winget install MiKTeX.MiKTeX)
# Font: TeX Gyre Termes (ships with MiKTeX, has Greek and maths - no DejaVu fallback in body text).

$ErrorActionPreference = 'Stop'

$src    = 'D:\BaiduSyncdisk\下载'
$work   = 'D:\ar_pilot\build'
$figs   = 'D:\ar_pilot\submission\figures'
$suppMd = 'D:\ar_pilot\supp_pandoc.md'
$pandoc = Join-Path $env:LOCALAPPDATA 'Pandoc\pandoc.exe'
$texbin = Join-Path $env:LOCALAPPDATA 'Programs\MiKTeX\miktex\bin\x64'
$env:PATH = "$texbin;$env:PATH"

New-Item -ItemType Directory -Force -Path "$work\figures" | Out-Null
Copy-Item (Join-Path $src 'manuscript_cited.md')  $work -Force
Copy-Item (Join-Path $src 'references_final.bib') $work -Force
Copy-Item (Join-Path $src 'biomed-central.csl') $work -Force -ErrorAction SilentlyContinue
foreach ($f in 'Figure1_expression','Figure2_calibration','Figure3_site_layer','Figure4_attrition') {
    Copy-Item (Join-Path $figs "$f.pdf") "$work\figures" -Force
}
foreach ($f in 'FigureS1_depth','FigureS2_nulls') {
    Copy-Item (Join-Path $figs "$f.png") "$work\figures" -Force
}

Push-Location $work

& $pandoc 'manuscript_cited.md' -o 'manuscript_final.pdf' --pdf-engine=xelatex --include-in-header='D://ar_pilot//build//bmc_header.tex' `
    --citeproc --bibliography='references_final.bib' --csl='biomed-central.csl' `
    -V mainfont="TeX Gyre Termes" -V fontsize=10pt `
    -V geometry:a4paper -V geometry:margin=2.54cm -V linestretch=2 `
    -H 'header.tex' --pdf-engine-opt=--enable-installer

Pop-Location

Push-Location 'D:\ar_pilot'

& $pandoc -f markdown-implicit_figures 'supp_pandoc.md' -o "$work\Supplementary_final.pdf" --pdf-engine=xelatex --include-in-header='D:/ar_pilot/build/supp_header.tex' `
    -V mainfont="TeX Gyre Termes" `
    -V geometry:a4paper -V geometry:margin=2cm -V linestretch=2 --columns=110 `
    -H 'D:\ar_pilot\build\supp9.tex' --pdf-engine-opt=--enable-installer

Pop-Location

foreach ($f in 'manuscript_final.pdf','Supplementary_final.pdf') {
    if (Test-Path "$work\$f") { Copy-Item "$work\$f" $src -Force; "OK  $f" } else { "FAILED  $f" }
}


# Editable manuscript file for journals that require DOCX (BMC accepts only
# DOC/DOCX/RTF/TeX as the manuscript file). Figures are embedded as PNG,
# because Word cannot display embedded PDF images.
$m = Get-Content "$work\manuscript_cited.md" -Raw -Encoding UTF8
$m = $m -replace 'figures/(Figure\d_\w+)\.pdf', 'figures/$1.png'
Set-Content -Path "$work\manuscript_docx.md" -Value $m -Encoding utf8
foreach ($f in 'Figure1_expression','Figure2_calibration','Figure3_site_layer','Figure4_attrition') {
    Copy-Item (Join-Path $figs "$f.png") "$work\figures" -Force
}
Push-Location $work
& $pandoc 'manuscript_docx.md' -o 'manuscript_final.docx' --citeproc --bibliography='references_final.bib' --csl='biomed-central.csl'
Pop-Location
if (Test-Path "$work\manuscript_final.docx") { Copy-Item "$work\manuscript_final.docx" $src -Force; "OK  manuscript_final.docx" } else { "FAILED  docx" }

