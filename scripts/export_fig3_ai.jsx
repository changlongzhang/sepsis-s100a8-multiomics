#target illustrator
(function(){
  app.userInteractionLevel=UserInteractionLevel.DONTDISPLAYALERTS;
  var b="D:/桌面/sepsis/文章/返修投稿新主图/";
  var d=app.open(new File(b+"Fig3.pdf"));
  var ao=new IllustratorSaveOptions();ao.pdfCompatible=true;ao.compressed=true;
  d.saveAs(new File(b+"Fig3.ai"),ao);
  d.close(SaveOptions.DONOTSAVECHANGES);
  var preview=new File(b+"Fig3_preview.png");if(preview.exists)preview.remove();
  app.userInteractionLevel=UserInteractionLevel.DISPLAYALERTS;
})();
