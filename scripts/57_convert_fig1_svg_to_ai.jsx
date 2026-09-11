#target illustrator

(function () {
    app.userInteractionLevel = UserInteractionLevel.DONTDISPLAYALERTS;
    var base = "C:/Users/ZCL/sepsis_workspace_link/review_revision/10_figure_sources/Fig1_S1_corrected_panels/";
    var doc = app.open(new File(base + "Fig1_CD_font_standardized_candidate.svg"));
    var options = new IllustratorSaveOptions();
    options.pdfCompatible = true;
    options.embedICCProfile = true;
    options.compressed = true;
    doc.saveAs(new File(base + "Fig1_CD_font_standardized_candidate.ai"), options);
    doc.close(SaveOptions.DONOTSAVECHANGES);
    app.userInteractionLevel = UserInteractionLevel.DISPLAYALERTS;
})();
