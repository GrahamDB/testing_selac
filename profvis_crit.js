
  // Calculate amount of time spent on each line of code. Returns nested objects
  // grouped by file, and then by line number.
  function getFileLineStats(prof, files) {
    // Drop entries with null or "" filename
    prof = prof.filter(function(row) {
      return row.filename !== null && row.filename !== "";
    });

    // Gather line-by-line file contents
    var fileLineStats = files.map(function(file) {
      // Create array of objects with info for each line of code.
      var lines = file.content.split("\n");
      var lineData = [];
      var filename = file.filename;
      var normpath = file.normpath;
      for (var i=0; i<lines.length; i++) {
        lineData[i] = {
          filename: filename,
          normpath: normpath,
          linenum: i + 1,
          content: lines[i],
          sumTime: 0,
          sumMem: 0,
          sumMemAlloc: 0,
          sumMemDealloc: 0
        };
      }

      return {
        filename: filename,
        lineData: lineData
      };
    });

    // Get timing data for each line
    var statsData = d3.nest()
      .key(function(d) { return d.filename; })
      .key(function(d) { return d.linenum; })
      .rollup(function(leaves) {
        var sumTime = leaves.reduce(function(sum, d) {
          // Add this node's time only if no ancestor node has the same
          // filename and linenum. This is to avoid double-counting times for
          // a line.
          var incTime = 0;
          if (!ancestorHasFilenameLinenum(d.filename, d.linenum, d.parent)) {
            incTime = d.endTime - d.startTime;
          }
          return sum + incTime;
        }, 0);

        var sumMem = leaves.reduce(function(sum, d) {
          return sum + d.sumMem;
        }, 0);

        var sumMemDealloc = leaves.reduce(function(sum, d) {
          return sum + d.sumMemDealloc;
        }, 0);

        var sumMemAlloc = leaves.reduce(function(sum, d) {
          return sum + d.sumMemAlloc;
        }, 0);

        return {
          filename: leaves[0].filename,
          linenum: leaves[0].linenum,
          sumTime: sumTime,
          sumMem: sumMem,
          sumMemAlloc: sumMemAlloc,
          sumMemDealloc: sumMemDealloc
        };
      })
      .entries(prof);

    // Insert the sumTimes into line content data
    statsData.forEach(function(fileInfo) {
      // Find item in fileTimes that matches the file of this fileInfo object
      var fileLineData = fileLineStats.filter(function(d) {
        return d.filename === fileInfo.key;
      })[0].lineData;

      fileInfo.values.forEach(function(lineInfo) {
        lineInfo = lineInfo.values;
        fileLineData[lineInfo.linenum - 1].sumTime = lineInfo.sumTime;
        fileLineData[lineInfo.linenum - 1].sumMem = lineInfo.sumMem;
        fileLineData[lineInfo.linenum - 1].sumMemDealloc = lineInfo.sumMemDealloc;
        fileLineData[lineInfo.linenum - 1].sumMemAlloc = lineInfo.sumMemAlloc;
      });
    });

    // Calculate proportional times, relative to the longest time in the data
    // set. Modifies data in place.
    var fileMaxTimes = fileLineStats.map(function(lines) {
      var lineTimes = lines.lineData.map(function(x) { return x.sumTime; });
      return d3.max(lineTimes);
    });

    var maxTime = d3.max(fileMaxTimes);

    fileLineStats.map(function(lines) {
      lines.lineData.map(function(line) {
        line.propTime = line.sumTime / maxTime;
      });
    });

    var totalMem = getTotalMemory(prof);

    fileLineStats.map(function(lines) {
      lines.lineData.map(function(line) {
        line.propMem = line.sumMem / totalMem;
        line.propMemDealloc = line.sumMemDealloc / totalMem;
        line.propMemAlloc = line.sumMemAlloc / totalMem;
      });
    });

    return fileLineStats;

    // Returns true if the given node or one of its ancestors has the given
    // filename and linenum; false otherwise.
    function ancestorHasFilenameLinenum(filename, linenum, node) {
      if (!node) {
        return false;
      }
      if (node.filename === filename && node.linenum === linenum) {
        return true;
      }
      return ancestorHasFilenameLinenum(filename, linenum, node.parent);
    }
  }

  function prepareProfData(prof, interval) {
    // Convert object-with-arrays format prof data to array-of-objects format
    var data = colToRows(prof);
    data = addParentChildLinks(data);
    data = consolidateRuns(data);
    data = applyInterval(data, interval);
    data = findCollapsedDepths(data);

    return data;
  }

  // Given the raw profiling data, convert `time` and `lastTime` fields to
  // `startTime` and `endTime`, and use the supplied interval. Modifies data
  // in place.
  function applyInterval(prof, interval) {
    prof.forEach(function(d) {
      d.startTime = interval * (d.time - 1);
      d.endTime = interval * (d.lastTime);
      delete d.time;
      delete d.lastTime;
    });

    return prof;
  }

  // Find the total time spanned in the data
  function getTotalTime(prof) {
    return d3.max(prof, function(d) { return d.endTime; }) -
           d3.min(prof, function(d) { return d.startTime; });
  }

  // Find the total memory spanned in the data
  function getTotalMemory(prof) {
    return d3.max(prof, function(d) { return d.memalloc; });
  }

  // Calculate the total amount of time spent in each function label
  function getAggregatedLabelTimes(prof) {
    var labelTimes = {};
    var tree = getProfTree(prof);
    calcLabelTimes(tree);

    return labelTimes;

    // Traverse the tree with the following strategy:
    // * Check if current label is used in an ancestor.
    //   * If yes, don't add to times for that label.
    //   * If no, do add to times for that label.
    // * Recurse into children.
    function calcLabelTimes(node) {
      var label = node.label;
      if (!ancestorHasLabel(label, node.parent)) {
        if (labelTimes[label] === undefined)
          labelTimes[label] = 0;

        labelTimes[label] += node.endTime - node.startTime;
      }

      node.children.forEach(calcLabelTimes);
    }

    // Returns true if the given node or one of its ancestors has the given
    // label; false otherwise.
    function ancestorHasLabel(label, node) {
      if (node) {
        if (node.label === label) {
          return true;
        }
        return ancestorHasLabel(label, node.parent);
      } else {
        return false;
      }
    }
  }


  // Given profiling data, add parent and child links to indicate stack
  // relationships.
  function addParentChildLinks(prof) {
    var data = d3.nest()
      .key(function(d) { return d.time; })
      .rollup(function(leaves) {
        leaves = leaves.sort(function(a, b) { return a.depth - b.depth; });

        leaves[0].parent = null;
        leaves[0].children = [];

        for (var i=1; i<leaves.length; i++) {
          leaves[i-1].children.push(leaves[i]);
          leaves[i].parent = leaves[i-1];
          leaves[i].children = [];
        }

        return leaves;
      })
      .map(prof);

    // Convert data from object of arrays to array of arrays
    data = d3.map(data).values();
    // Flatten data
    return d3.merge(data);
  }


  // Given profiling data, consolidate consecutive blocks for a flamegraph.
  // This function also assigns correct parent-child relationships to form a
  // tree of data objects, with a hidden root node at depth 0.
  function consolidateRuns(prof) {
    // Create a special top-level leaf whose only purpose is to point to its
    // children, the items at depth 1.
    var topLeaf = {
      depth: 0,
      parent: null,
      children: prof.filter(function(d) { return d.depth === 1; })
    };

    var tree = consolidateTree(topLeaf);
    var data = treeToArray(tree);
    // Remove the root node from the flattened data
    data = data.filter(function(d) { return d.depth !== 0; });
    return data;

    function consolidateTree(tree) {
      var leaves = tree.children;
      leaves = leaves.sort(function(a, b) { return a.time - b.time; });

      // Collapse consecutive leaves, with some conditions
      var startLeaf = null;  // leaf starting this run
      var lastLeaf = null;   // The last leaf we've looked at
      var newLeaves = [];
      var collectedChildren = [];
      var sumMem = 0;
      var sumMemDealloc = 0;
      var sumMemAlloc = 0;

      // This takes the start leaf, end leaf, and the set of children for the
      // new leaf, and creates a new leaf which copies all its properties from
      // the startLeaf, except lastTime and children.
      function addNewLeaf(startLeaf, endLeaf, newLeafChildren, sumMem, sumMemDealloc, sumMemAlloc) {
        var newLeaf = $.extend({}, startLeaf);
        newLeaf.lastTime = endLeaf.time;
        newLeaf.parent = tree;
        newLeaf.children = newLeafChildren;

        // Recurse into children
        newLeaf = consolidateTree(newLeaf);

        // Aggregate memory from this consolidation batch and their children
        aggregateMemory(newLeaf, sumMem, sumMemDealloc, sumMemAlloc);

        newLeaves.push(newLeaf);
      }

      function aggregateMemory(leaf, sumMem, sumMemDealloc, sumMemAlloc) {
        leaf.sumMem = sumMem;
        leaf.sumMemDealloc = sumMemDealloc;
        leaf.sumMemAlloc = sumMemAlloc;
        if (leaf.children) {
          leaf.children.forEach(function(child) {
            leaf.sumMem += child.sumMem ? child.sumMem : 0;
            leaf.sumMemDealloc += child.sumMemDealloc ? child.sumMemDealloc : 0;
            leaf.sumMemAlloc += child.sumMemAlloc ? child.sumMemAlloc : 0;
          });
        }
      }

      for (var i=0; i<leaves.length; i++) {
        var leaf = leaves[i];

        if (i === 0) {
          startLeaf = leaf;
          sumMem = sumMemAlloc = sumMemDealloc = 0;
        } else if (leaf.label !== startLeaf.label ||
                   leaf.filename !== startLeaf.filename ||
                   leaf.linenum !== startLeaf.linenum ||
                   leaf.depth !== startLeaf.depth)
        {
          addNewLeaf(startLeaf, lastLeaf, collectedChildren, sumMem, sumMemDealloc, sumMemAlloc);

          collectedChildren = [];
          startLeaf = leaf;
          sumMem = sumMemAlloc = sumMemDealloc = 0;
        }

        sumMem += leaf.meminc;
        sumMemDealloc += Math.min(leaf.meminc, 0);
        sumMemAlloc += Math.max(leaf.meminc, 0);
        collectedChildren = collectedChildren.concat(leaf.children);
        lastLeaf = leaf;
      }

      // Add the last one, if there were any at all
      if (i !== 0) {
        addNewLeaf(startLeaf, lastLeaf, collectedChildren, sumMem, sumMemDealloc, sumMemAlloc);
      }

      tree.children = newLeaves;
      return tree;
    }

    // Given a tree, pull out all the leaves and put them in a flat array
    function treeToArray(tree) {
      var allLeaves = [];

      function pushLeaves(leaf) {
        allLeaves.push(leaf);
        leaf.children.forEach(pushLeaves);
      }

      pushLeaves(tree);
      return allLeaves;
    }
  }


  // Given profiling data with parent-child information, get the root node.
  function getProfTree(prof) {
    if (prof.length === 0)
      return null;

    // Climb up to the top of the tree
    var node = prof[0];
    while (node.parent) {
      node = node.parent;
    }
    return node;
  }
