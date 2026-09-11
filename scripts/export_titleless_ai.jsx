#target illustrator
(function(){
  app.userInteractionLevel=UserInteractionLevel.DONTDISPLAYALERTS;
  var b="D:/桌面/sepsis/文章/返修投稿新主图/";
  var names=["Fig1","Fig4","Fig5","Fig6","Fig7","Fig8","S2_fig","S3_fig","S4_fig","S5_fig","S6_fig","S7_fig","S8_fig"];
  for(var i=0;i<names.length;i++){
    var n=names[i], d=app.open(new File(b+n+".pdf"));
    var arial=null;
    try{arial=app.textFonts.getByName("ArialMT");}catch(e1){try{arial=app.textFonts.getByName("Arial");}catch(e2){}}
    if(arial!==null){for(var j=0;j<d.textFrames.length;j++)d.textFrames[j].textRange.characterAttributes.textFont=arial;}
    var ao=new IllustratorSaveOptions();ao.pdfCompatible=false;ao.compressed=true;
    d.saveAs(new File(b+n+".ai"),ao);
    d.close(SaveOptions.DONOTSAVECHANGES);
  }
  app.userInteractionLevel=UserInteractionLevel.DISPLAYALERTS;
})();
