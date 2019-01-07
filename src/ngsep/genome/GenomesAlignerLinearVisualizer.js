var margin = {top: 45, right: 80, bottom: 30, left: 80},
    width = 1300 - margin.left - margin.right,
    height = 900 - margin.top - margin.bottom;
  
var dimensions = [
  {
    name: "geneStart",
    scale: d3.scale.linear().range([height, 0]),
    type: "number"
  },
  {
    name: "geneEnd",
    scale: d3.scale.linear().range([height, 0]),
    type: "number"
  },
  {
    name: "chromosome",
    scale: d3.scale.ordinal().rangePoints([0, height]),
    type: "string"
  },
  {
    name: "geneId",
    scale: d3.scale.ordinal().rangePoints([0, height]),
    type: "string"
  },
  {
    name: "geneIdG2",
    scale: d3.scale.ordinal().rangePoints([0, height]),
    type: "string"
  },
  {
    name: "chromosomeG2",
    scale: d3.scale.ordinal().rangePoints([0, height]),
    type: "string"
  },
  {
    name: "geneStartG2",
    scale: d3.scale.linear().range([height, 0]),
    type: "number"
  },
  {
    name: "geneEndG2",
    scale: d3.scale.linear().range([height, 0]),
    type: "number"
  },
];
  
var dimensionss = [
  {
    name: "averageA",
    scale: d3.scale.linear().range([height, 0]),
    type: "number"
  },
  {
    name: "averageB",
    scale: d3.scale.linear().range([height, 0]),
    type: "number"
  }
];

var x = d3.scale.ordinal().rangePoints([0, width]),
    y = {},
    dragging = {};

var c20 = d3.scale.category20();

var line = d3.svg.line(),
    axis = d3.svg.axis(),
    background,
    foreground;

var svg = d3.select("body").append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
        .call(d3.behavior.zoom().on("zoom", function() {
      svg.attr("transform", "translate("+d3.event.translate+")"+ " scale("+ d3.event.scale+")")
    }))
    // .call(d3.behavior.zoom().on("zoom", function() {
      // svg.attr("transform", "scale("+ d3.event.scale+")")
   //  }))
  .append("g")
    .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

d3.tsv("InputFileLCS.tsv", function(error, lcs) 
{
  d3.tsv("InputFileGenome1.tsv", function(error, crmsA)
  {
    d3.tsv("InputFileGenome2.tsv", function(error, crmsB)
    {
      crmsA[0].startt = 0
      crmsA[0].endd = crmsA[0].Length
      for(var i = 1; i<crmsA.length; i++){
        crmsA[i].startt = crmsA[i-1].endd
        crmsA[i].endd = ""+(Number(crmsA[i].Length)+Number(crmsA[i].startt))
      }
      chrmsA = toObject(crmsA);
      
      crmsB[0].startt = 0
      crmsB[0].endd = crmsB[0].Length
      for(var i = 1; i<crmsB.length; i++){
        crmsB[i].startt = crmsB[i-1].endd
        crmsB[i].endd = ""+(Number(crmsB[i].Length)+Number(crmsB[i].startt))
      }
      chrmsB = toObject(crmsB);
      
      for(var j = 0; j<lcs.length; j++){
        genomeA = lcs[j].chromosome;
        start = Number(lcs[j].geneStart);
        end = Number(lcs[j].geneEnd);
        lcs[j].averageA = Number(chrmsA[genomeA].startt)+(start+end)/2;
    
        genomeB = lcs[j].chromosomeG2;
        startG2 = Number(lcs[j].geneStartG2);
        endG2 = Number(lcs[j].geneEndG2);
        lcs[j].averageB = Number(chrmsB[genomeB].startt)+(startG2+endG2)/2;
    
      }
      // Extract the list of dimensions and create a scale for each.
      x.domain(dimensionss.map(function(d) { return d.name; }))
      dimensionss.forEach(function(dimension) {
      dimension.scale.domain(dimension.type === "number"
        ? d3.extent(lcs, function(d) { return +d[dimension.name]; })
        : lcs.map(function(d) { return d[dimension.name]; }).sort());
      y[dimension.name] = dimension.scale;
  });

  // Add grey background lines for context.
  //background = svg.append("g")
  //    .attr("class", "background")
  //  .selectAll("path")
  //    .data(lcs)
  //  .enter().append("path")
  //    .attr("d", path);

  // Add blue foreground lines for focus.
  //foreground = svg.append("g")
  //    .attr("class", "foreground")
  //  .selectAll("path")
  //    .data(lcs)
  //  .enter().append("path")
  //    .attr("d", path)
  //  	This function assigns color to genes aligning to the same chromosome 
  //    .attr("stroke",function(d){return c20(d.chromosome)});
      

 //Show LCS only
d3.select("#LCSbutton").on("click", showLCS);  
function showLCS() {
  // Add colored foreground lines for focus.
  foreground = svg.append("g")
      .attr("class", "foreground")
    .selectAll("path")
      .data(lcs)
    .enter().append("path")
      .attr("d", path)
  	  .attr("stroke", function(d) {if(d.type=="L") return c20(d.chromosome)})
}
 //Show multiple hits
d3.select("#Multiplebutton").on("click", showMultiple);
function showMultiple() {
  // Add colored foreground lines for focus.
  foreground = svg.append("g")
      .attr("class", "foreground")
    .selectAll("path")
      .data(lcs)
    .enter().append("path")
      .attr("d", path)
  	  .attr("stroke", function(d) {if(d.type=="M") return c20(d.chromosome)})
}
  
//Show unique hits
d3.select("#Uniquesbutton").on("click", showUniques);
function showUniques() {
  // Add colored foreground lines for focus.
  foreground = svg.append("g")
      .attr("class", "foreground")
    .selectAll("path")
      .data(lcs)
    .enter().append("path")
      .attr("d", path)
  	  .attr("stroke", function(d) {if(d.type=="U") return c20(d.chromosome)})
}      
      

  // Add a group element for each dimension.
  var g = svg.selectAll(".dimension")
      .data(dimensionss)
    .enter().append("g")
      .attr("class", "dimension")
      .attr("transform", function(d) {  return "translate(" + x(d.name) + ")"; })
      .call(d3.behavior.drag()
        .origin(function(d) { return {x: x(d.name)}; })
        .on("dragstart", function(d) {
          dragging[d.name] = x(d.name);
          background.attr("visibility", "hidden");
        })
        .on("drag", function(d) {
          dragging[d.name] = Math.min(width, Math.max(0, d3.event.x));
          foreground.attr("d", path);
          dimensionss.sort(function(a, b) { return position(a) - position(b); });
          x.domain(dimensionss.map(function(d) { return d.name; }));
          g.attr("transform", function(d) { return "translate(" + position(d) + ")"; })
        })
        .on("dragend", function(d) {
          delete dragging[d.name];
          transition(d3.select(this)).attr("transform", "translate(" + x(d.name) + ")");
          transition(foreground).attr("d", path);
          background
              .attr("d", path)
            .transition()
              .delay(500)
              .duration(0)
              .attr("visibility", null);
        }));

 tempA = []
    tempA[0]=0
    teempA = {0:""}
    for(var i = 1; i<=crmsA.length; i++){
        tempA[i]=parseInt(crmsA[i-1].endd)
        teempA[crmsA[i-1].endd]=crmsA[i-1].Name
    }
    console.log(tempA)
    tempB = []
    tempB[0]=0
    teempB = {0:""}
    for(var i = 1; i<=crmsB.length; i++){
        tempB[i]=parseInt(crmsB[i-1].endd)
        teempB[crmsB[i-1].endd]=crmsB[i-1].Name
    }
  // Add an axis and title.
  g.append("g")
      .attr("class", "axis")
      .each(function(d) {
          if(d.name == "averageA"){
            d3.select(this)
              .call(axis.scale(y[d.name]).orient("left").tickValues(tempA).tickFormat(function(d){return teempA[d]}))
              .attr("class","axis averageA")
              .selectAll(".averageA .tick text")
              .on('click', function (d, i) {
                console.log(teempA[d])
              })
              .style("text-anchor", "start")
              .attr('transform', function(d){
                if(d!=0) {
                  var yy = y["averageA"]((parseInt(chrmsA[teempA[d]].startt)+parseInt(d))/2)
                  var xx = y["averageA"](parseInt(d))
                  return "translate(-65,"+(yy-xx)+")"
                }
                return "translate(0,10)"
              }) 
          } else {
            d3.select(this)
              .call(axis.scale(y[d.name]).orient("right").tickValues(tempB).tickFormat(function(d){return teempB[d]}))
              .attr("class","axis averageB")
              .selectAll(".averageB .tick text")
              .on('click', function (d, i) {
                console.log(d)
              })
              .style("text-anchor", "start")
              .attr('transform', function(d){
                if(d!=0) {
                  var yy = y["averageB"]((parseInt(chrmsB[teempB[d]].startt)+parseInt(d))/2)
                  var xx = y["averageB"](parseInt(d))
                  return "translate(0,"+(yy-xx)+")"
                }
                return "translate(0,10)"
              })
              // .selectAll(".tick text")
              // .style("text-anchor", "start")
              // .attr("x", axisTextAdjust);  
          }
        })
    .append("text")
      .style("text-anchor", "right")
      .attr("y", -20)
     // .text("Genome")
      //.text(function(d) { return d.name; });

  function adjustTextLabels(selection) {
    selection.selectAll('.major text')
      ;
  }



  // Add and store a brush for each axis.
  g.append("g")
      .attr("class", "brush")
      .each(function(d) {
        d3.select(this).call(y[d.name].brush = d3.svg.brush().y(y[d.name]).on("brushstart", brushstart).on("brush", brush));
      })
    .selectAll("rect")
      .attr("x", -8)
      .attr("width", 16);
  })})
});

function position(d) {
  var v = dragging[d.name];
  return v == null ? x(d.name) : v;
}

function transition(g) {
  return g.transition().duration(500);
}

// Returns the path for a given data point.
function path(d) {
  //console.log(d)
  //var a = line(dimensions.map(function(p) { 
    //console.log(p)
    //return [position(p), y[p](d[p])]; }));
  //console.log(dimensions.map(function(p) { return [position(p), y[p](d[p])]; }))
  //return a;
  linee = line(dimensionss.map(function(dimension) {
    return [x(dimension.name), dimension.scale(d[dimension.name])];
  }));
  return linee;
}

function brushstart() {
  d3.event.sourceEvent.stopPropagation();
}

// Handles a brush event, toggling the display of foreground lines.
function brush() {
  var actives = dimensionss.filter(function(p) {
    return !y[p.name].brush.empty(); }),
      extents = actives.map(function(p) { return y[p.name].brush.extent(); });
  foreground.style("display", function(d) {
    return actives.map(a => a.name).every(function(p, i) {
      var w = isNaN(parseInt(d[p]))?y[p](d[p]):d[p]
      return extents[i][0] <= w && w <= extents[i][1];
    }) ? null : "none";
  });
}
  
function toObject(arr) {
  var rv = {};
  for (var i = 0; i < arr.length; ++i)
    rv[arr[i].Name] = arr[i];
  return rv;
}
