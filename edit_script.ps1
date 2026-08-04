$path = "R\plot.R"
$content = [System.IO.File]::ReadAllText($path)

# Edit 1: change paste0(tag, "\n", pat) to paste(tag, pat)
$old1 = '      paste0(tag, "\n", pat)'
$new1 = '      paste(tag, pat)  # Single line: "R PAT_001" — works with 45° strip'
if ($content.Contains($old1)) {
    $content = $content.Replace($old1, $new1)
    Write-Output "Edit 1 applied"
} else {
    Write-Output "Edit 1 not found (already applied?)"
}

# Edit 2: wrap drug in brackets
$old2 = '  y_lab <- if (!is.null(drug)) paste0("Predicted Viability (z-score)\n", drug)' + "`n" + '           else "Predicted Viability (z-score)"'
$new2 = '  y_lab <- if (!is.null(drug)) paste0("Predicted Viability (z-score)\n[", drug, "]")' + "`n" + '           else "Predicted Viability (z-score)"'
if ($content.Contains($old2)) {
    $content = $content.Replace($old2, $new2)
    Write-Output "Edit 2 applied"
} else {
    Write-Output "Edit 2 not found (already applied?)"
}

# Edit 3: add y-axis title margin
$old3 = '      axis.title = element_text(size = rel(0.95)),'
$new3 = '      axis.title = element_text(size = rel(0.95)),' + "`n" + '      axis.title.y = element_text(margin = margin(r = 8, unit = "pt")),'
if ($content.Contains($old3)) {
    $content = $content.Replace($old3, $new3)
    Write-Output "Edit 3 applied"
} else {
    Write-Output "Edit 3 not found (already applied?)"
}

[System.IO.File]::WriteAllText($path, $content)
Write-Output "Done."
