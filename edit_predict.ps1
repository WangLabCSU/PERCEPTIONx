$predictPath = 'c:\Users\Lenovo\Desktop\PERCEPTION\R\predict.R'
$content = Get-Content -Path $predictPath -Raw

$oldText = "  # Subset expression matrix to model features`n  dataset_FOI <- data.frame(`n    dataset[feature_match, ],`n    row.names = rownames(dataset)[feature_match]`n  )`n  dataset_FOI <- data.frame(t(dataset_FOI))`n  predict(model, dataset_FOI)`n}"

$newText = "  # Subset expression matrix to model features`n  orig_colnames <- colnames(dataset)[feature_match]`n  dataset_FOI <- data.frame(`n    dataset[feature_match, ],`n    row.names = rownames(dataset)[feature_match]`n  )`n  dataset_FOI <- data.frame(t(dataset_FOI))`n  # data.frame() mangles special chars (e.g. ""@@"" -> "".."") via make.names().`n  # Restore original colnames so predictions carry the proper clone keys.`n  colnames(dataset_FOI) <- orig_colnames`n  predict(model, dataset_FOI)`n}"

if ($content.Contains($oldText)) {
    $content = $content.Replace($oldText, $newText)
    Set-Content -Path $predictPath -Value $content -NoNewline
    Write-Output "predict.R: SUCCESS"
} else {
    Write-Output "predict.R: FAILED - old text not found"
}
