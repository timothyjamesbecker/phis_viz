//globals------------------------------------------------------------------
var categories = null,
    time_range = [new Date('01/01/1920'), new Date('01/01/2020')],
    full_tree_width_factor = 0.95,
    less_tree_width_factor = 1.5,
    tree_width_factor = 0.95,
    color_range = [180,0],
    rad_mult = 50,
    rad_off = 14,
    trunks = 4,
    power_factor = 0.5;

function get_categories(pt){
    var C = {},
        cs = [];
    for(var i = 0; i<pt.length; i++) {
        var k = parseInt(pt[i].Category);
        if (k in C) {
            C[k] += 1;
        } else {
            C[k] = 1;
            cs.push(parseInt(pt[i].Category));
        }
    }
    var min = 10000000,
        max = -10000000;
    for(var i = 0; i < cs.length; i++){
        if(cs[i]<min){ min = cs[i]; }
        if(cs[i]>max){ max = cs[i]; }
    }
    return max-min+1;
}

function normalize_samples(st){
    var min = 1000000,
        max = 0,
        diff = 0;
    for(var i = 0; i < st.length; i++){
        st[i].RelAbund = parseFloat(st[i].RelAbund);
        if(st[i].RelAbund > max){
            max = st[i].RelAbund;
        }
        if(st[i].RelAbund < min){
            min = st[i].RelAbund;
        }
    }

    console.log(max);
    console.log(min);
    if(max-min>0) {
        diff = max - min;
        for (var i = 0; i < st.length; i++) {
            st[i].RelAbund = Math.pow((st[i].RelAbund - min)/diff,power_factor);
        }
    }
    return st;
}

function category_to_dist(d,m){
    if('category' in d){
        var c = parseInt(d.category);
        d.category = [];
        for(var i = 0; i < m; i++){
            if(i==c){ d.category.push(1.0); }
            else{     d.category.push(0.0); }
        }
    }
}

//given dist [3.0,4.0,0.0],[h1,h2] => h
function dist_to_color(dist,H){
    var sum = d3.sum(dist),
        d = dist.length,
        a = 0.0;
    for(var i = 0; i < d; i++){
        dist[i] = dist[i]/sum;
        a = a + (i/(d-1)*dist[i]);     //right side 1.0]
    }
    if(H[1]>H[0]){
        return a*(H[1]-H[0])+H[0];
    }else {
        return (1.0-a)*(H[0]-H[1])+H[1];
    }
}

function hsla_color_gradient(dist,H,alpha){
    return 'hsla(' + dist_to_color(dist,H) + ',100%,50%,'+alpha.toString()+')';
}

function get_temporal_set(tree){
    var temporal_set = {},
        temporal_list = [];
    function dft(T){
        if(!('children' in T || '_children' in T)){//leaf node
            if('samples' in T.data) {
                for (var j = 0; j < T.data.samples.length; j++) {
                    var temporal = T.data.samples[j].Collection_Date.toString();
                    if (temporal in temporal_set) {
                        temporal_set[temporal] = temporal_set[temporal] + 1;
                    } else {
                        temporal_set[temporal] = 1;
                    }
                }
            }
        }else{
            if('children' in T){
                for(var i = 0; i < T.children.length; i++){
                    dft(T.children[i]);
                }
            }else{
                for(var i = 0; i < T._children.length; i++){
                    dft(T._children[i]);
                }
            }
        }
    }
    dft(tree);
    for(var temporal in temporal_set){
        var date = new Date(temporal);
        var dates = [date.getMonth()+1,date.getDate(),date.getFullYear()];
        temporal_list.push(dates.join('/'));
    }
    return temporal_list.sort(function(a,b){ return (new Date(a)) > (new Date(b)); });
}

//st is the sample table array of objects that has every data point
//s_idx is a species index that maps the names to the alph position
//sids are the unique ids you will encounter in the st
function sample_table_temporal_analysis(st,s_idx){
    var S = {},
        stats = {'min':1.0,'max':0.0,'sum':0.0,'n':st.length,
                 'd_min':1.0,'d_max':0.0,'d_sum':0.0, 'center':new_vector(s_idx),
                 'diffs':{},'norm_diffs':{},'sites':{},'t_start':time_range[1],'t_stop':time_range[0]};

    function new_vector(idx){
        var m = [];
        for(var i in idx){ m.push(0.0); }
        return m;
    }

    //[1] populate the S and stats
    for(var i = 0; i < st.length; i++){
        var l = st[i].SID,
            date = st[i].Collection_Date,
            t = [date.getMonth()+1,date.getDate(),date.getFullYear()].join('/'),
            s = st[i].Taxon_name,
            v = st[i].RelAbund;
        var t_time = new Date(t);
        if(t_time<stats.t_start){ stats.t_start = t_time; }
        if(t_time>stats.t_stop){  stats.t_stop  = t_time; }
        if(v<stats.min){ stats.min = v; } //global min
        if(v>stats.max){ stats.max = v; } //global max
        stats.sum += v;                   //global sum

        if(l in S){
            if(!(t in S[l])){ S[l][t] = new_vector(s_idx); } //new time row
            S[l][t][s_idx[s]] = v;  //set species value
            stats.sites[l].sum += v; //update site stats
            stats.sites[l].n += 1.0;
            if(v<stats.sites[l].min){ stats.sites[l].min = v; }
            if(v>stats.sites[l].max){ stats.sites[l].max = v; }
        }else{ //new site, start a new time/species matrix
            S[l] = {};                    //new location
            S[l][t] = new_vector(s_idx);  //new time row
            S[l][t][s_idx[s]] = v;        //set the species value
            stats.sites[l] = {'min':v,'max':v,'sum':v,'n':1.0,
                              'd_min':1.0,'d_max':0.0,'d_sum':0.0,
                              'center': new_vector(s_idx)}; //init stats
            stats.diffs[l]      = {};
            stats.norm_diffs[l] = {};
        }
    }

    //[2] center and diffs of global and each site using any/all temporal slices
    for(var l in S){
        var n = 1.0*Object.keys(S[l]).length,
            m = 1.0*Object.keys(s_idx).length;
        for(var t in S[l]){
            for(var s in s_idx){
                stats.sites[l].center[s_idx[s]] += S[l][t][s_idx[s]]/n;
            }
        }
        for(var t in S[l]){
            stats.diffs[l][t] = 0.0;
            for(var s in s_idx){
                var x = (S[l][t][s_idx[s]]-stats.sites[l].center[s_idx[s]])/m;
                if(x<0.0){ x *= -1.0; }
                stats.diffs[l][t] += x;
            }
            stats.sites[l].d_sum += stats.diffs[l][t];
            if(stats.diffs[l][t] < stats.sites[l].d_min){ stats.sites[l].d_min = stats.diffs[l][t]; }
            if(stats.diffs[l][t] > stats.sites[l].d_max){ stats.sites[l].d_max = stats.diffs[l][t]; }
        }
        stats.d_sum += stats.sites[l].d_sum;
        if(n>1){
            if(stats.sites[l].d_min < stats.d_min){ stats.d_min = stats.sites[l].d_min; }
            if(stats.sites[l].d_max > stats.d_max){ stats.d_max = stats.sites[l].d_max; }
        }
    }

    //[3] normalize differences now that you have analyzed all the data
    var x = 0;
    for(var l in S){
        for(var t in S[l]){
            if((stats.sites[l].d_max-stats.sites[l].d_min)>0.0){
                stats.norm_diffs[l][t] = (stats.diffs[l][t]-stats.d_min)/(stats.d_max-stats.d_min);
                if(stats.norm_diffs[l][t]<0.0){ stats.norm_diffs[l][t] = 0.0; }
                if(stats.norm_diffs[l][t]>1.0){ stats.norm_diffs[l][t] = 1.0; }
            }
            x += 1;
        }
    }
    stats.S = S;
    return stats;
}

//change to site coordinate center
function site_center(G){
    var xy = [0.0,0.0],
        n = 0;
    for(var i = 0; i < G.features.length; i++){
        if(G.features[i].geometry.type=='Point') {
            xy[1] = xy[1] + G.features[i].geometry.coordinates[0];
            xy[0] = xy[0] + G.features[i].geometry.coordinates[1];
            n = n + 1.0;
        }
    }
    if(n >= 1.0) {
        xy[0] = xy[0] / n;
        xy[1] = xy[1] / n;
    }
    //console.log(xy);
    return xy;
}

function distance_vector(sid,time,S,method){
    var D = {};
    if(method=='bray-curtis'){
        var a_i = S[sid][time],
            a_sum = 0;
        for(var i in a_i){ a_sum += a_i[i]; } //total sum of site
        for(var l in S){ //now look at all other sites
            D[l] = 0.0; //will take only the maximum similarity for sites
            for(var t in S[l]){
                if(t!=time) { //remove the actual sample, but this leaves room for other time samples at that site
                    var b_sum = 0.0,
                        d_sum = 0.0,
                        d_ij = 0.0;
                    for (var j in S[l][t]) {
                        b_sum += S[l][t][j];
                        d_sum += Math.abs(a_i[j]-S[l][t][j]);
                    }
                    if (a_sum + b_sum > 0.0) {
                        d_ij = 1.0 - d_sum / (a_sum + b_sum);
                        if (d_ij > D[l]) { D[l] = d_ij; } //keep only the max
                    }
                }
            }
        }
    }
    return D;
}

function time_distance_vector(sid_source,time_source,s_id,S,method){
    var D = {};
    if(method=='bray-curtis'){
        var a_i = S[sid_source][time_source],
            a_sum = 0.0;
        for(var i in a_i){ a_sum += a_i[i]; } //total sum of source site
        for(var t in S[s_id]){ //now look at the target sites time slots
            D[t] = 0.0; //will take only the maximum similarity for sites
            var b_sum = 0.0,
                d_sum = 0.0;
            for(var j in S[s_id][t]){
                b_sum += S[s_id][t][j];
                d_sum += Math.abs(a_i[j]-S[s_id][t][j]);
            }
            if(a_sum+b_sum>0.0){
                D[t] = 1.0-d_sum/(a_sum+b_sum);
            }
        }
    }
    return D;
}

//given the phis tree, return the leaf nodes as a sorted list
function get_species_map_from_tree(tree){
    var species = {},
        species_idx = {};
    function dft(tree){
        if(!('children' in tree || '_children' in tree)){
            if(!(tree.name in species)) {
                species[tree.name] = 1;
            }else{
                species[tree.name] += 1;
            }
        }else{
            if('children' in tree){
                for(var i = 0; i < tree.children.length; i++){
                    dft(tree.children[i]);
                }
            }else{
                for(var i = 0; i < tree._children.length; i++){
                    dft(tree._children[i]);
                }
            }
        }
    }
    dft(tree);
    var order = Object.keys(species).sort();
    for(var s in species) {
        for (var i = 0; i < order.length; i++) {
            if(s==order[i]) {
                species_idx[s] = i;
                break;
            }
        }
    }
    return species_idx;
}

function get_species_data(species,tree){
    var data = {};
    function dft(T){
        if(!('children' in T || '_children' in T)){
            if(T.name == species) {
               data = T.data;
            }
        }else{
            if('children' in T){
                for(var i = 0; i < T.children.length; i++){
                    dft(T.children[i]);
                }
            }else{
                for(var i = 0; i < T._children.length; i++){
                    dft(T._children[i]);
                }
            }
        }
    }
    dft(tree);
    return data;
}

function get_species_map_from_pt(pt){
    var species = {},
        species_idx = {};
    for(var i = 0; i < pt.length; i++){
        if(pt[i].Taxon_name in species){
            species[pt[i].Taxon_name] = 1;
        }else{
            species[pt[i].Taxon_name] += 1;
        }
    }
    var order = Object.keys(species).sort();
    for(var s in species) {
        for (var i = 0; i < order.length; i++) {
            if(s==order[i]) {
                species_idx[s] = i;
                break;
            }
        }
    }
    return species_idx;
}

/*take a diatom tree table and calculate the viz tree
  D is a parsed d3.csv data object:
  [{Phylum:p1,Class:c1,Order:o1,Family:f1,Genus:g1,Taxon_name:t1,Short_name:s1,Category:l1},...]
if trim is set, find the straddle point where the tree starts to split and go one back
*/
function phylo_table_to_viz_tree(D,trim){
    var S = {},
        T = {},
        R = [];
    for(var i = 0; i < D.length; i++){
        var r = [];
        Object.keys(D[i]).forEach(function(d){
            r.push(D[i][d]);
            if(D[i][d] in S){ S[D[i][d]]++; }
            else{             S[D[i][d]]=0; }
        });
        R.push(r);
    }

    //make d an object so you can build it our accross several data sets
    function insert(V,row,i,d){
        if(i>=row.length-3){ //leaf node => path=[leaf]
            var t = {'name':row[i],'data':d};
            category_to_dist(t.data,categories);
            t.data.descendants=0; //update freq by definition of leaf
            t.data.short_name=row[row.length-2];
            V.children.push(t); //no children on the leaf nodes
        }
        else{ //inner nodes
            if(i<=0){ //root node
                if(!('name' in V)){
                    V.name = row[i];
                    V.data = {'descendants':1};
                    V.children = [];
                }else{
                    V.data.descendants++; //update freq
                }
                insert(V,row,i+1,d);
            }else {
                var v = V.children.map(function(x){ return x.name; }).indexOf(row[i]);
                if(v<0) { //wasn't found among the children
                    var t = {'name': row[i], 'data':{'descendants':1},'children':[]};
                    V.children.push(t); //parent inserts the link to new inner node
                    insert(t,row,i+1,d);
                }else{
                    V.children[v].data.descendants++;  //update freq
                    insert(V.children[v],row,i+1,d);
                }
            }
        }
    }
    //read in the rows of data and update the viz tree structure
    for(var r = 0; r < R.length; r++){
        insert(T,R[r],0,{'category':R[r][R[r].length-1]});
    }
    T.meta = {'n':T.data.descendants};
    if(trim){
        function diverge(T){
            if(!('children' in T)){ //leaf node
                return T; //return just the leaf
            }else{
                if(T.children.length>1){
                    return T;
                }else{
                    return diverge(T.children[0]);
                }
            }
        }
        T = diverge(T);
    }
    return T;
}

/*given a Taxon_name and the sample_table, return all unique sids*/
function phlyo_to_site_ids(st,name,normalize){
    var sids = {},
        min = 1.0,
        max = 0.0;
    for(var i = 0; i < st.length; i++){
        if(st[i].Taxon_name==name){
            if(st[i].SID in sids){
                sids[st[i].SID][0] += st[i].RelAbund;
                sids[st[i].SID][1] += 1;
            }else{
                sids[st[i].SID] = [st[i].RelAbund,1];
            }
        }
    }
    for(var sid in sids){
        sids[sid][0] /= 1.0*sids[sid][1]; // average
        if(sids[sid][0]>=max){ max = sids[sid][0]; }
        if(sids[sid][0]<=min){ min = sids[sid][0]; }
    }
    for(var sid in sids){ //normalize [0.0,1.0]
        if(max-min>0.0 && normalize) {
            sids[sid][0] = (sids[sid][0]-min)/(max - min);
            if (sids[sid][0] < 0.0) {
                sids[sid][0] = 0.0;
            }
            if (sids[sid][0] > 1.0) {
                sids[sid][0] = 1.0;
            }
        }
    }
    return sids;
}

function sample_table_join_viz_tree(S,T,tr,sid,prop){
    var time_set
    function insert(T,name,d){
        if(T.name==name){ //should be a leaf
            if('samples' in T.data){ T.data.samples.push(d); }
            else{ T.data.samples = [d]; }
            T.data.samples.sort(function(a,b){ return a.Collection_Date > b.Collection_Date; }); //ascending by Date(s)
        }else{
            if('children' in T) {
                for (var i = 0; i < T.children.length; i++) {
                    insert(T.children[i],name,d);
                }
            }
        }
    }

    for(var i = 0; i < S.length; i ++){
        S[i].SID = parseInt(S[i].SID);
        S[i].RelAbund = parseFloat(S[i].RelAbund);
        S[i].Collection_Date = new Date(S[i].Collection_Date);
        if(S[i].SID==sid && (S[i].Collection_Date>=tr[0] && S[i].Collection_Date<=tr[1])) {
            insert(T,S[i].Taxon_name,S[i]);
        }
    }

    if(prop){
        function propagation(V,c){
            if(!('children' in V)) { //leaf node
                if('samples' in V.data){
                    V.data.RelAbund = 1e-12;
                    for(var j = 0; j < V.data.samples.length;j++){
                        V.data.RelAbund = V.data.RelAbund+V.data.samples[0].RelAbund;
                    }
                    V.data.RelAbund = V.data.RelAbund/V.data.samples.length;
                } else{
                    V.data.RelAbund = 1e-12;
                }
                return {'RelAbund':V.data.RelAbund,'category':V.data.category};
            }else{         //root or inner nodes
                V.data.RelAbund = 0.0;
                V.data.category = [];
                for(var j = 0; j < c; j++){ V.data.category.push(0.0); }
                for(var i = 0; i < V.children.length; i++){
                    var data = propagation(V.children[i],c);
                    V.data.RelAbund = V.data.RelAbund+data.RelAbund;
                    for(var j = 0; j < c; j++) {
                        V.data.category[j] = V.data.category[j]+data.RelAbund*data.category[j];
                    }
                }
                return {'RelAbund':V.data.RelAbund,'category':V.data.category};
            }
        }
        propagation(T,categories); //3 categories => array of 3 => [1.0,0.0,0.0]
    }
    //insert some meta data on a sorted time set to graph.
    return T;
}

function phlyo_tree_map_graph(phylo_data_url,sample_data_url,phylo_id,
                              site_geojson_url,map_id,mapbox) {
    var attach_id = document.getElementById('phis_viz');
    var margin = {top:   attach_id.clientHeight/45,
                  bottom:attach_id.clientHeight/45,
                  right: attach_id.clientWidth/45,
                  left:  attach_id.clientWidth/45};
    var width = attach_id.clientWidth/tree_width_factor - (margin.left + margin.right);
    var height = attach_id.clientHeight - (margin.top + margin.bottom);
    var i = 0,
        duration = 300,
        root = {'x0':height/2,'y0':0,'sim':false},
        selected_phis_node = null,
        last_phis_color = null,
        temporal_list = [],
        selected_time_slice = null,
        time_state = {pt:null,st:null,tl:null,sid:null,
                      mode:'distance',s_idx:null,refresh:true},
        hover_slice = null,
        last_hover_slice = null,
        last_time_slice = null,
        phylo_zoom = false,
        phylo_collapse = false,
        phylo_table = [],
        sample_table = [],
        sample_stats = {},
        s_idx = null,
        geojsonFeature = null,
        lastsitesLayer = null,
        lastlinesLayer = null,
        selected_geo_json = {'last':null},
        mymap = null;

    var diagonal = d3.svg.diagonal()
        .projection(function (d) {
            return [d.y, d.x];
        });

    var svg = d3.select('#'+phylo_id).append("svg")
        .attr("width", width + margin.right + margin.left)
        .attr("height", height + margin.top + margin.bottom)
        .call(d3.behavior.zoom().on("zoom",function(){
            if(phylo_zoom) {
                var xy = d3.event.translate;
                xy[0] = xy[0] + 2 * margin.left;
                xy[1] = xy[1] + 2 * margin.top;
                svg.attr("transform", "translate(" + xy[0] + "," + xy[1] + ")" + " scale(" + d3.event.scale + ")")
            }
        }))
        .append("g")
        .attr("transform", "translate(" + 2*margin.left + "," + margin.top + ")");

    /*
        dt is the diaton tree structure, st is the site table data structure to insert into the dt
        tr is the time_range object to query in the st and the sid is the sid to query in the st
     */
    function refresh_tree(pt,st,tr,sid){
        var x0 = root.x0;
        var y0  = root.y0;
        var sim = root.sim;
        root = phylo_table_to_viz_tree(pt,true);
        root = sample_table_join_viz_tree(st,root,tr,sid,true); //could select right here...
        root.x0 = x0;
        root.y0 = y0;
        root.sim = sim;
        function collapse(d) {
            if (d.children) {
                d._children = d.children;
                d._children.forEach(collapse);
                d.children = null;
            }
        }
        //root.children.forEach(collapse); //collapse the children
        update(root);
    }

    function refresh_temporal_controls_timeline(){
        d3.select('#temporal_bar').html('');

        var pt   = time_state.pt,
            st   = time_state.st,
            tl   = time_state.tl,
            sid  = time_state.sid,
            mode = time_state.mode;
        if(tl!=null) {
            function update_time_line() {
                var D = null,
                    sim = selected_geo_json['sim'], //set to null when nothing is selected
                    node = d3.select(selected_phis_node).node();
                if (root.sim && node != null) { //you have a selected_phis_node and simularity mode is on...
                    var parent_node = d3.select(selected_phis_node).node().parentNode;
                    var parent = d3.select(parent_node).select('circle');
                    var species = $(parent_node).attr('name');
                    // console.log(species);

                    var time = tl[tl.length-1];
                    var target = sid;
                    if(selected_geo_json['time']!=null){
                        time = d3.select(selected_geo_json['time'][0][0]).attr('date');
                    }
                    if(sim!=null){
                        target = sim.target.feature.properties.SID;
                    }
                    D = time_distance_vector(target,time,sid,sample_stats.S,'bray-curtis');
                } else {
                    if(node!=null) { //you have species mode here, must pull the data from sample_table?
                        var parent_node = d3.select(selected_phis_node).node().parentNode;
                        var parent = d3.select(parent_node).select('circle');
                        var species = $(parent_node).attr('name');
                        var s_idx = get_species_map_from_tree(root);
                        D = {};
                        for(var i = 0; i < tl.length; i++){
                            D[tl[i]] = sample_stats.S[sid][tl[i]][s_idx[species]];
                        }
                        var data = get_species_data(species,root);
                        // console.log('trying to locate the classification vector to then derive the correct color');
                        // console.log(data);
                        if(Object.keys(data).length>0) {
                            // console.log(hsla_color_gradient(data.category,color_range,0.75));
                            last_phis_color = hsla_color_gradient(data.category, color_range, 0.75);
                        }
                    }else{
                        if (selected_time_slice == null) {
                            D = time_distance_vector(sid,tl[tl.length-1],sid,sample_stats.S,'bray-curtis');
                        } else {
                            var time = d3.select(selected_time_slice[0][0]).attr('date');
                            if(!(time in sample_stats.S[sid])){
                                D = time_distance_vector(sid,tl[tl.length-1],sid,sample_stats.S,'bray-curtis');
                            }else {
                                D = time_distance_vector(sid, time, sid, sample_stats.S, 'bray-curtis');
                            }
                        }
                    }
                }

                var time_points = [];
                for (var i = 0; i < tl.length; i++) { //gather up your circle data
                    time_points.push({'date': new Date(tl[i]), 'value': D[Object.keys(D)[i]]});
                }

                t_svg.selectAll('circle').data([]).exit().remove();
                var timeline = t_svg.selectAll('circle')
                    .data(time_points)
                    .enter()
                    .append('circle')
                    .attr('r', function (d,i) { //will be the radius
                        if (selected_time_slice == null && i == time_points.length - 1) {
                            selected_time_slice = d3.select(this);
                            selected_geo_json['time'] = selected_time_slice;
                        }
                        if(isNaN(d.value)){ return rad_off;  }
                        else{ return rad_off + rad_mult * d.value; }
                    })
                    .attr('cy', function (d) {
                        return h / 2 - 2;
                    })
                    .attr('cx', function (d) {
                        return x(d.date);
                    })
                    .attr('date', function (d) {
                        return [d.date.getMonth() + 1, d.date.getDate(), d.date.getFullYear()].join('/');
                    })
                    .style('fill', function (d) {
                        var c = 'rgba(255,255,255,0.5)';
                        var t = [d.date.getMonth()+1,d.date.getDate(),d.date.getFullYear()].join('/');
                        if (t == d3.select(selected_time_slice)[0][0].attr('date')) {
                            c = 'rgba(255,255,255,0.9)';
                        }
                        if(!root.sim && node!=null){ //species mode is on, we need to get the color
                            c = last_phis_color      //not robust against temporal absence now
                            if(c.startsWith('rgba')) {
                                rgba = c.split('rgba(')[1].split(')')[0].split(',')
                                rgba[3] = 0.5;
                                if (t == d3.select(selected_time_slice)[0][0].attr('date')) {
                                    rgba[3] = 0.9;
                                }
                                c = 'rgba(' + rgba.join(',') + ')';
                            }else{
                                hsla = c.split('hsla(')[1].split(')')[0].split(',')
                                hsla[3] = 0.5;
                                if (t == d3.select(selected_time_slice)[0][0].attr('date')) {
                                    hsla[3] = 0.9;
                                }
                                c = 'hsla(' + hsla.join(',') + ')';
                            }
                        }
                        return c;
                    })
                    .on('mouseover', function (d) { //change color and display time point
                        var time = [d.date.getMonth() + 1, d.date.getDate(), d.date.getFullYear()].join('/');
                        var nav_text = $('#nav_site').text().split(' > ')[0].split(' | ')[0];
                        var nav = $('#nav_site').text(nav_text+' | TID='+time);
                        d3.select(this).style({'stroke': 'rgba(255,255,255,0.5)', 'stroke-width': '1vh'});
                    })
                    .on('mouseout', function (d) { //change back to unselected color and turn off display time point
                        var nav = $('#nav_site').text($('#nav_site').text().split(' | TID=')[0]);
                        d3.select(this).style({'stroke': 'rgba(255,255,255,0.0)', 'stroke-width': '1vh'});
                        d3.select(this).text('');
                    })
                    .on("click", function (d) { //update the selection pointers and refresh the tree
                        selected_time_slice = d3.select(this);
                        time_state.refresh = true;
                        refresh_tree(pt, st, [d.date, d.date], sid);
                        last_time_slice = selected_time_slice;
                        update_time_line();
                    });
            }

            var margin = {
                left: 24,
                right: 24,
                top: 0,
                bottom: 0
            };
            var w = attach_id.clientWidth - (margin.left + margin.right);
            var h = attach_id.clientHeight * 0.065;

            var x = d3.time.scale()
                .range([0, w - (margin.left + margin.right)])
                .domain(d3.extent([{'date': sample_stats.t_start}, {'date': sample_stats.t_stop}],
                    function (d) {
                        return d.date;
                    }));

            var xAxis = d3.svg.axis()
                .scale(x)
                .orient('top')
                .ticks(7);

            var t_svg = d3.select('#temporal_bar').append('svg')
                .attr('width', w)
                .attr('height', h)
                .style({'position': 'fixed', 'bottom': '2vh'})
                .append('g')
                .attr('transform', 'translate(' + margin.left + ',' + margin.top + ')');

            // Add the X Axis
            t_svg.append("g")
                .attr("class", "x axis")
                .attr("transform", 'translate(0,' + h + ')')
                .style({
                    'color': 'rgba(255,255,255,0.9)', 'fill': 'none',
                    'stroke': 'rgba(255,255,255,0.9)', 'stroke-width': 1,
                    'font-family': 'Lato', 'font-size': '2vh'
                })
                .call(xAxis);
            //set the last_time_slice to the newest
            update_time_line();
            if(time_state.refresh) {
                refresh_tree(pt, st, [new Date(tl[tl.length-1]), new Date(tl[tl.length-1])], sid);
                time_state.refresh = false;
            }
            last_time_slice = selected_time_slice;
        }
    }

    $.getJSON(site_geojson_url,function(geo_data) {
        d3.csv(phylo_data_url, function (phylo_error, phylo_data) {
            d3.csv(sample_data_url, function (sample_error, sample_data) {
                if (phylo_error) throw phylo_error;
                if (sample_error) throw sample_error;
                phylo_table = phylo_data;
                categories = get_categories(phylo_table);
                sample_table = normalize_samples(sample_data);

                //map time-----------------------------------------
                geojsonFeature = geo_data; //load your geojson data
                var attach_id = document.getElementById('phis_viz');
                var width = 3 * (attach_id.clientWidth / 2.0);
                var height = attach_id.clientHeight;
                $(map_id).css({'height': height, 'width': width});

                mymap = new L.Map(map_id);
                mymap.setView(site_center(geojsonFeature), 9)
                mymap.doubleClickZoom.disable();
                //basemap loading via OSM or mapbox
                if (mapbox) {
                    L.tileLayer('https://api.tiles.mapbox.com/v4/{id}/{z}/{x}/{y}.png?access_token=pk.eyJ1IjoibWFwYm94IiwiYSI6ImNpejY4NXVycTA2emYycXBndHRqcmZ3N3gifQ.rJcFIG214AriISLbB6B5aw', {
                        maxZoom: 16,
                        minZoom: 9,
                        attribution: 'Map data &copy; <a href="http://openstreetmap.org">OpenStreetMap</a> contributors, ' +
                        '<a href="http://creativecommons.org/licenses/by-sa/2.0/">CC-BY-SA</a>, ' +
                        'Imagery © <a href="http://mapbox.com">Mapbox</a>',
                        id: 'mapbox.streets'
                    }).addTo(mymap);
                } else {
                    // create the tile layer with correct attribution
                    //var osmUrl = 'http://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png';
                    var osmUrl = 'https://server.arcgisonline.com/ArcGIS/rest/services/Canvas/World_Light_Gray_Base/MapServer/tile/{z}/{y}/{x}'
                    var osmAttrib = 'Map data © <a href="https://openstreetmap.org">OpenStreetMap</a> contributors';
                    var osm = new L.TileLayer(osmUrl, {minZoom: 9, maxZoom: 16, attribution: osmAttrib});
                    mymap.addLayer(osm);
                }

                var selected_options = {
                    radius: rad_mult/6.0,
                    fillColor: "rgb(250,250,250)",
                    color: "rgb(255,255,255)",
                    opacity: 1.0,
                    fillOpacity: 0.5
                };

                var unselected_options = {
                    radius: rad_mult/8.0,
                    fillColor: "rgb(25,25,25)",
                    color: "rgb(200,200,200)",
                    opacity: 0.6,
                    fillOpacity: 0.75
                };

                var hover_options = {
                    radius: rad_mult/12.0,
                    fillColor: "rgb(75,75,75)",
                    color: "rgb(255,255,255)",
                    opacity: 1.0,
                    fillOpacity: 0.75
                };

                function onEachFeature(feature, layer) {
                    layer.on('click',function(e){
                        var site = e.target.feature.properties;
                        if(selected_geo_json['last']==null){ //no sites was on, start a new select
                            e.target.setStyle(selected_options);
                        }else{ //a site was last selected
                            var last = selected_geo_json['last'].target.feature.properties
                            if(site.SID!=last.SID){ //new site selected
                                // console.log('new site selected');
                                selected_geo_json['last'].target.setStyle(unselected_options);
                                e.target.setStyle(selected_options);
                                //a new site was selected, see if simulation if on
                                var sim = selected_geo_json['sim'], //set to null when nothing is selected
                                    node = d3.select(selected_phis_node).node();
                                if(sim==null){ //reset the selected_time_slice
                                    console.log('resetting time slice');
                                    selected_time_slice = null;
                                }else{
                                    console.log('new site, sim mode on!')
                                }
                            }else{ //already have it selected, turn off tree
                                // console.log('already was selected')
                                if(e.target.options.fillColor==selected_options.fillColor) {
                                    e.target.setStyle(unselected_options);
                                    site = {'SID':-1,'Station_Name':''};
                                }else{
                                    e.target.setStyle(selected_options);
                                }
                            }
                        }
                        selected_geo_json['last'] = e;
                        time_state.refresh = true;
                        refresh_tree(phylo_table,sample_table,time_range,site.SID);
                        if(site.SID>=0) {
                            var nav = $('#nav_site').text('SID=' + site.SID + ' : ' + site.Station_Name.toUpperCase());
                        }else{
                            var nav = $('#nav_site').text('phis_viz');
                        }
                        temporal_list = get_temporal_set(root);
                        $('#temporal_bar').html('');
                        if(temporal_list.length>0) { //make a set of date buttons that will invoke samples
                            time_state.pt = phylo_table;
                            time_state.st = sample_table;
                            time_state.tl = temporal_list;
                            time_state.sid = site.SID;
                            time_state.mode = 'distance';
                            refresh_temporal_controls_timeline();
                        }
                    });

                    layer.on('mouseover',function(e){
                        var site = e.target.feature.properties;
                        var s_idx = get_species_map_from_tree(root);
                        var extra_metric = '';
                        var sim = '';

                        // console.log(e.target);
                        // console.log(sample_stats.S[site.SID]);
                        // console.log(selected_phis_node);
                        // console.log(s_idx);

                        if(selected_phis_node!=null){ //make sure you have leaf or root picked first
                            if(!root.sim) { //species mode use relative abundance
                                // console.log('species mode used...');
                                var node = d3.select(selected_phis_node).node().parentNode;
                                var species = d3.select(node).attr('name');
                                if (true) {//check to see if you are in sim mode or species
                                    var sum = 0.0,
                                        n = 0;
                                    for (var t in sample_stats.S[site.SID]) {
                                        sum += sample_stats.S[site.SID][t][s_idx[species]];
                                        n += 1;
                                    }
                                    extra_metric = ' ' + d3.format('.5f')(sum / n) + ' mean relative abundance';
                                }
                            }
                        }

                        var last = {'SID':-1};
                        var nav_text = $('#nav_site').text();
                        if(selected_geo_json['last']!=null) {
                            last = selected_geo_json['last'].target.feature.properties;
                        }
                        if(site.SID!=last.SID) {
                            e.target.setStyle(hover_options);
                            var nav = $('#nav_site').text(nav_text+' > '+site.Station_Name.toUpperCase()+extra_metric);
                        }
                    });

                    layer.on('mouseout',function(e){
                        var site = e.target.feature.properties;
                        var last = {'SID':-1};
                        var nav_text = $('#nav_site').text();
                        if(selected_geo_json['last']!=null) {
                            last = selected_geo_json['last'].target.feature.properties;
                        }
                        if(site.SID!=last.SID){
                            e.target.setStyle(unselected_options);
                            var nav_text = $('#nav_site').text($('#nav_site').text().split(' > ')[0]);
                        }
                    });
                }
                L.geoJSON(geojsonFeature,{
                    pointToLayer: function (feature, latlng) {
                        return L.circleMarker(latlng, unselected_options);
                    },
                    onEachFeature: onEachFeature
                }).addTo(mymap);

                lastsitesLayer = L.geoJSON({type: "FeatureCollection",features:[]},{
                    pointToLayer: function (feature, latlng) {
                        return L.circleMarker(latlng, unselected_options);
                    },
                    onEachFeature: onEachFeature
                }).addTo(mymap);

                mymap.on('zoomend', function() {
                    var zoom_value = mymap.getZoom();
                });
                time_state.refresh = true;
                refresh_tree(phylo_table,sample_table,time_range,-1); //init is a 0.0 relabund blank tree
                //load and attach the per site temporal analysis once and then access it later
                s_idx = get_species_map_from_pt(phylo_table);
                sample_stats = sample_table_temporal_analysis(sample_table,s_idx);
            });
        });
    });

    $('#map_toggle_button').click(function(){
        var w = parseInt($('#map_side').css('width').split('px')[0]);
        var v = document.getElementById('phis_viz').clientWidth;
        if(w<=0.1){
            // console.log('map is open...');
            tree_width_factor = full_tree_width_factor;
            $('#map_side').css('width','49%');
        }
        if(w/v > 0.45 && w/v < 0.55){
            // console.log('map is closed...');
            tree_width_factor = less_tree_width_factor;
            $('#map_side').css('width','0%');
        }
        update(root);
    });

    $('#data_check_button').click(function(){
        //[1] do a file check on the data
        var phlyo_report = {},
            sample_report = {},
            geo_report = {};
        d3.csv(phylo_data_url, function (phylo_error, phylo_data) {
            d3.csv(sample_data_url, function (sample_error, sample_data) {
                d3.json(site_geojson_url,function(geo_error, geo_data) {
                    if(phylo_error){ phylo_report.error = phylo_error; }
                    if(sample_error){ sample_report.error = sample_error; }
                    if(geo_error){ geo_report.error = geo_error; }
                });
            });
        });
        //[2] launch an alert with the report for each data section
    });



    function update(source) {
        var margin = {top:   attach_id.clientHeight/45,
                      bottom:attach_id.clientHeight/45,
                      right: attach_id.clientWidth/45,
                      left:  attach_id.clientWidth/45};
        var w = parseInt($('#map_side').css('width').split('px')[0]);
        if(w<=1){
            // console.log('something weird happend');
            var width = attach_id.clientWidth - (margin.left + margin.right);
        }else {
            var width = attach_id.clientWidth/tree_width_factor - (margin.left + margin.right);
        }

        var height = attach_id.clientHeight - (margin.top + margin.bottom);
        var tree = d3.layout.tree().size([height,width+margin.left+margin.right]);
        d3.select(phylo_id).select('svg')
            .attr("width", width + margin.right + margin.left)
            .attr("height", height + margin.top + margin.bottom);

        // Compute the new tree layout.
        var nodes = tree.nodes(root).reverse(),
            links = tree.links(nodes);
        // Normalize for fixed-depth.
        nodes.forEach(function (d) {
            d.y = d.depth * width/8;
        });
        // Update the nodes…
        var node = svg.selectAll("g.node")
            .data(nodes, function (d) {
                return d.id || (d.id = ++i);
            });
        // Enter any new nodes at the parent's previous position.
        var nodeEnter = node.enter().append("g")
            .attr("class", "node")
            .attr("transform", function(d){
                if(source==null || !('y0' in source) || !('x0' in source)){ /*console.log(source);*/ }
                return "translate(" + source.y0 + "," + source.x0 + ")";
            })
            .attr('name',function(d){ return d.name });
        nodeEnter.append("circle")
            .attr("r", 1e-6)
            .on('mouseover',mouse_over)
            .on('mouseout',mouse_out)
            .on("click", left_click)        //toggle children off and on
        //.on('contextmenu', d3.contextMenu(menu));  //context menu or add node?

        nodeEnter.append("text")
            .attr("x", function (d) {
                return d.children || d._children ? 0 : 10;
            })
            .attr("y", function (d) {
                return d.children || d._children ? -1*(d.data.RelAbund*rad_mult+rad_off) : 0;
            })
            .attr("dy", ".35em")
            .style('font-size','0.75vw')
            .attr("text-anchor", function (d) {
                return d.children || d._children ? "end" : "start";
            })
            .text(function (d) {
                if(d.data.descendants<1){ return d.data.short_name; }
                var stub = d.name.split(' ');
                if(stub.length>1) { stub = stub.splice(1).join(' '); }
                else{ stub = stub[0]; }
                return d.children ? '' : stub;
            })
            .style("fill-opacity", 1e-6);

        // Transition nodes to their new position.
        var nodeUpdate = node.transition()
            .duration(duration)
            .attr("transform", function (d) {
                if(!('y' in d) || !('x' in d)){ /*console.log(d);*/ }
                return "translate(" + d.y + "," + d.x + ")";
            });

        nodeUpdate.select("circle")
            .attr("r", function(d){
                return d._children ? rad_mult*d.data.RelAbund+rad_off:(rad_mult-rad_off)*d.data.RelAbund+rad_off;
            })
            .style("fill", function(d) {
                var c = hsla_color_gradient(d.data.category,color_range,0.75);
                if (d.data.RelAbund < 1e-9){
                    c = 'rgba(255,255,255,0.5)';
                }
                if(selected_phis_node!=null){
                    var selected = d3.select(d3.select(selected_phis_node).node().parentNode);
                    var selected_name = selected.attr('name');
                    if(d.name==selected_name) {
                        c = 'rgba(120,120,120,0.5)';
                    }
                }
                return d._children ? "rgba(0,0,0,0.5)":c;
            })

        nodeUpdate.select("text")
            .attr("x", function (d) {
                return d.children || d._children ? 0 : rad_off;
            })
            .attr("y", function (d) {
                return d.children || d._children ? -1*(d.data.RelAbund*rad_mult+rad_off) : 0;
            })
            .attr("dy", ".35em")
            .attr("text-anchor", function (d) {
                return d.children || d._children ? "end" : "start";
            })
            .style("fill-opacity", 1);

        // Transition exiting nodes to the parent's new position.
        var nodeExit = node.exit().transition()
            .duration(duration)
            .attr("transform", function (d) {
                if(source==null || !('y' in source) || !('x' in source)){ }//console.log(source); }
                else {
                    return "translate(" + source.y + "," + source.x + ")";
                }
            })
            .remove();

        nodeExit.select("circle")
            .attr("r", 1e-6);
        nodeExit.select("text")
            .attr("x", function (d) {
                return d.children || d._children ? 0 : rad_off;
            })
            .attr("y", function (d) {
                return d.children || d._children ? -1*(d.data.RelAbund*rad_mult+rad_off) : 0;
            })
            .attr("dy", ".35em")
            .attr("text-anchor", function (d) {
                return d.children || d._children ? "end" : "start";
            })
            .style("fill-opacity", 1e-6);

        // Update the links…
        var link = svg.selectAll("path.link")
            .data(links, function (d) {
                return d.target.id;
            });

        // Enter any new links at the parent's previous position.
        link.enter().insert("path", "g")
            .attr("class", "link")
            .attr("d", function (d) {
                var o = {x: source.x0, y: source.y0};
                return diagonal({source: o, target: o});
            });
        // mouseover for edges/paths can be truned off or on
        // .on('mouseover',function(d){
        //     d3.select(this).style("stroke-width",function(d){
        //         return d.target.data.RelAbund*150.0+8;
        //     })
        //     d3.select('#nav_data').html(d3.format(".5f")(d.target.data.RelAbund));
        // })
        // .on('mouseout',function(d){
        //     d3.select(this).style("stroke-width",function(d){
        //         return d.target.data.RelAbund*100.0+4;
        //     })
        //     d3.select('#nav_data').html(d3.format(".5f")(''));
        // });

        // Transition links to their new position.
        link.transition()
            .duration(duration)
            .attr("d", diagonal)
            .style("stroke-width",function(d){
                return d.target.data.RelAbund*trunks*rad_mult+rad_off;
            })
            .style("stroke",function(d){
                if(d.target.data.RelAbund<1e-9){
                    var c = 'rgba(0,0,0,0.0);'
                }else{
                    //console.log(dist_to_color(d.target.data.category,color_range))
                    var c = 'hsla('+dist_to_color(d.target.data.category,color_range)+',100%,50%,0.75)';
                }
                return c;
            });

        // Transition exiting nodes to the parent's new position.
        link.exit().transition()
            .duration(duration)
            .attr("d", function (d) {
                var o = {x: source.x, y: source.y};
                return diagonal({source: o, target: o});
            })
            .remove();

        // Stash the old positions for transition.
        nodes.forEach(function (d) {
            d.x0 = d.x;
            d.y0 = d.y;
        })
    }

    function mouse_over(d){
        var self = d3.select(this);
        self.attr("r", function(d){
            //console.log(d);
            return (rad_mult+rad_off)*d.data.RelAbund+rad_off;
        })
            .attr("y", function (d) {
                return d.children || d._children ? -1*(d.data.descendants*(rad_mult-rad_off)) : 0;
            });
        var display = ' relative abundance for '+d.name;
        d3.select('#nav_data').html(d3.format(".5f")(d.data.RelAbund)+display);
    }

    function mouse_out(d){
        var self = d3.select(this);
        self.attr("r", function (d) {
            return (rad_mult-rad_off)*d.data.RelAbund + rad_off;
        });
        // if('children' in d || '_children' in d) {
        //     d3.select(this.parentNode).select('text')
        //         .text(function (d) {
        //             return '';
        //         })
        // }
        d3.select('#nav_data').html('');
    }

    function left_click(d){
        if(phylo_collapse){ //needs more work
            if (d.children) {
                d._children = d.children;
                d.children = null;
            }else{
                d.children = d._children;
                d._children = null;
            }
        }

        var presence_options = {
            radius: rad_mult/4.0,
            fillColor: "rgb(255,255,255)",
            color: "rgb(50,50,50)",
            opacity: 0.2,
            fillOpacity: 0.6
        };
        var similarity_options = {
            radius: rad_mult/4.0,
            fillColor: "rgb(255,255,255)",
            color: "rgb(50,50,50)",
            opacity: 0.8,
            fillOpacity: 0.2
        };
        var sites = phlyo_to_site_ids(sample_table,d.name,false);
        var sitesFeature = {type: "FeatureCollection",features:[]};

        if(d.data.descendants<1){//leaf click...
            // console.log('sim is off');
            root.sim = false; //turn off sim mode when a leaf clicked
            for(var i = 0; i < geojsonFeature.features.length; i++){
                if(geojsonFeature.features[i].properties.SID in sites){
                    sitesFeature.features.push(geojsonFeature.features[i]);
                }
            }
            mymap.removeLayer(lastsitesLayer);
            //keep track of the last node
            var node = d3.select(this).node();
            var parent_node = null;
            var species = '';
            if(selected_phis_node!=null){
                var parent_node = d3.select(selected_phis_node).node().parentNode;
                species = $(parent_node).attr('name');
            }
            // console.log(species);
            // console.log(d);

            //issue when you have species mode on and then select a new site and then toggle off species mode
            if(d.name==species) { //toggle off species mode
                // console.log('species mode off');
                var parent_node = d3.select(this).node().parentNode;
                var parent = d3.select(parent_node).select('circle');
                parent.style('fill',last_phis_color);

                lastsitesLayer = L.geoJSON({type: "FeatureCollection",features:[]},{
                    pointToLayer: function (feature, latlng) {
                        return L.circleMarker(latlng, presence_options);
                    }//,onEachFeature: onEachFeature
                }).addTo(mymap);
                lastsitesLayer.bringToBack();
                selected_phis_node = null;
                update(d);
            }else{ //turn on the species mode layer
                // console.log('species mode on');
                if(selected_phis_node!=null) {
                    var parent_node = d3.select(selected_phis_node).node().parentNode;
                    var parent = d3.select(parent_node).select('circle');
                    parent.style('fill', last_phis_color);
                }
                selected_phis_node = node;
                parent_node = d3.select(this).node().parentNode;
                parent = d3.select(parent_node).select('circle');
                last_phis_color = parent.style('fill');
                lastsitesLayer = L.geoJSON(sitesFeature,{
                    pointToLayer: function (feature, latlng) {
                        presence_options.radius = rad_off+0.5*rad_mult*sites[feature.properties.SID][0];
                        presence_options.fillColor = 'hsla(' + dist_to_color(d.data.category, color_range) + ',100%,50%,0.75)';
                        return L.circleMarker(latlng, presence_options);
                    }//,onEachFeature: onEachFeature
                }).addTo(mymap);
                lastsitesLayer.bringToBack();
                update(d);
            }
            refresh_temporal_controls_timeline();
        }else{ //inner nodes
            if(lastlinesLayer!=null){//toggle off the selected geo stuff
                mymap.removeLayer(lastlinesLayer);
                lastlinesLayer = null;
            }

            if(d.parent==null){ //root node click
                root.sim = !root.sim; //complete the toggle and now check if its on

                if(selected_phis_node!=null) { //reset the old selected leaf color
                    var parent_node = d3.select(selected_phis_node).node().parentNode;
                    var parent = d3.select(parent_node).select('circle');
                    parent.style('fill', last_phis_color);
                }
                parent_node = d3.select(this).node().parentNode;
                parent = d3.select(parent_node).select('circle');
                last_phis_color = parent.style('fill');
                parent.style('fill', 'rgba(120,120,120,0.5)');

                //turn off the site layer
                mymap.removeLayer(lastsitesLayer);
                lastsitesLayer = L.geoJSON({type: "FeatureCollection",features:[]},{
                    pointToLayer: function (feature, latlng) {
                        return L.circleMarker(latlng, presence_options);
                    }//,onEachFeature: onEachFeature
                }).addTo(mymap);
                lastsitesLayer.bringToBack();
                refresh_temporal_controls_timeline();
                if(root.sim && selected_geo_json!=null && selected_geo_json['last']!=null){ //turn on the layer!
                    // console.log('sim is on and site is selected!')

                    //if the selected_geo_json['last']!=slected_geo_json['sim'] => update last_time_slice?
                    //if the last_time_slice isn't in sample_stats.S[last.SID] => reset to new site_id
                    selected_geo_json['sim'] = selected_geo_json['last'];
                    selected_geo_json['time'] = last_time_slice;

                    var node = d3.select(this).node();
                    if(selected_phis_node!=null) { //reset the old selected leaf color
                        var parent_node = d3.select(selected_phis_node).node().parentNode;
                        var parent = d3.select(parent_node).select('circle');
                        parent.style('fill', last_phis_color);
                    }

                    selected_phis_node =  d3.select(this).node();
                    if(selected_geo_json['last']!=null) { //have a SID selected, now get the time slice, etc...
                        var last = selected_geo_json['last'].target.feature.properties;
                        var time = d3.select(last_time_slice[0][0]).attr('date');

                        if(!(time in sample_stats.S[last.SID])) { //forgot to click a new timepoint
                            var temp_list = Object.keys(sample_stats.S[last.SID]);
                            selected_geo_json['time'] = null;
                            time = temp_list[temp_list.length-1];
                            console.log('selecting time automatically');

                            selected_time_slice = null;
                            last_time_slice = null;
                            time_state.pt = phylo_table;
                            time_state.st = sample_table;
                            time_state.tl = temp_list;
                            time_state.sid = last.SID;
                            time_state.mode = 'distance';
                            time_state.refresh = true;
                            refresh_temporal_controls_timeline();
                        }

                        //dump text into into the id=nav_data span
                        $('#nav_search').text('SID='+last.SID+' | TID='+time);

                        var D = distance_vector(last.SID,time,sample_stats.S,'bray-curtis',0.75);
                        for(var i = 0; i < geojsonFeature.features.length; i++){
                            if(geojsonFeature.features[i].properties.SID in D){
                                sitesFeature.features.push(geojsonFeature.features[i]);
                            }
                        }

                        var filter = 0.75,
                            lastPoint = null,
                            simFeature = [],
                            simPoints  = [];
                        for(var i in sitesFeature.features){
                            var feature = sitesFeature.features[i];
                            if(feature.properties.SID==last.SID){
                                lastPoint = [feature.geometry.coordinates[1],feature.geometry.coordinates[0]];
                            }else{
                                if(D[feature.properties.SID]>=filter){
                                    // console.log({0:feature.properties.SID,1:D[feature.properties.SID]});
                                    simFeature.push(feature);
                                    simPoints.push([feature.geometry.coordinates[1],feature.geometry.coordinates[0]]);
                                }
                            }
                        }
                        var lineStyle = {color: 'rgba(0,0,0,0.1)'};
                        var lineLayers = [];
                        for(var i in simPoints){
                            var line = L.polyline([lastPoint,simPoints[i]],lineStyle);
                            line.bringToBack();
                            lineLayers.push(line);
                        }
                        lastlinesLayer = L.featureGroup(lineLayers).addTo(mymap);
                        lastlinesLayer.bringToBack();

                        function onEachSim(feature, layer) {
                            layer.on('mouseover', function (e) {
                                var site = e.target.feature.properties;
                                var nav_search = $('#nav_search').text();
                                var nav = $('#nav_search').text(d3.format('1.5f')(site.sim)+' max similarity to '+nav_search);
                            });
                            layer.on('mouseout',function(e){
                                var nav_search = $('#nav_search').text().split(' similarity to ');
                                var nav = $('#nav_search').text(nav_search[nav_search.length-1]);
                            });
                        }

                        lastsitesLayer = L.geoJSON(sitesFeature,{
                            pointToLayer: function (feature, latlng) {
                                feature.properties.sim = D[feature.properties.SID];
                                if(D[feature.properties.SID]>=filter) {
                                    similarity_options.radius = rad_off+0.5*rad_mult*D[feature.properties.SID];
                                }else{
                                    similarity_options.radius = 0;
                                }
                                if(feature.properties.SID==last.SID) {
                                    // similarity_options.fillColor = 'rgb(255,0,0)';
                                    // similarity_options.color = 'rgb(255,0,0)';
                                    similarity_options.fillOpacity = 1.0;
                                    similarity_options.opacity = 1.0;
                                    similarity_options.radius = rad_off+0.5*rad_mult* D[feature.properties.SID];
                                }else{
                                    similarity_options.fillColor = 'rgb(100,100,100)';
                                    similarity_options.color = 'rgb(50,50,50)';
                                    similarity_options.fillOpacity = 0.2;
                                    similarity_options.opacity = 0.6;
                                }
                                return L.circleMarker(latlng, similarity_options);
                            },
                            onEachFeature: onEachSim
                        }).addTo(mymap);
                        lastsitesLayer.bringToBack();
                    }
                    update(d);
                }else{ //d.sim is false so you need to turn it off
                    $('#nav_search').text('');
                    // console.log('sim is off or a site was not selected');
                    if(selected_geo_json!=null && selected_geo_json['last']!=null){//turn off the layer
                        var sid = selected_geo_json['last'].target.feature.properties.SID;
                        var temp_list = Object.keys(sample_stats.S[sid]);

                        selected_geo_json['sim'] = null;
                        selected_geo_json['time'] = null;
                        selected_phis_node = null;

                        lastsitesLayer = L.geoJSON(sitesFeature,{
                            pointToLayer: function (feature, latlng) {
                                similarity_options.radius = rad_off+0.5*rad_mult*D[feature.properties.SID];
                                return L.circleMarker(latlng, similarity_options);
                            }//,onEachFeature: onEachFeature
                        }).addTo(mymap);
                        lastsitesLayer.bringToBack();

                        var time = d3.select(last_time_slice[0][0]).attr('date');

                        // console.log(temp_list);
                        // console.log(time);

                        var t = false;
                        for(var i in temp_list){
                            if(temp_list[i]==time){
                                t = true;
                                break;
                            }
                        }
                        // console.log(t);

                        if(!t) { //forgot to click a new timepoint
                            selected_geo_json['time'] = null;
                            time = temp_list[temp_list.length-1];
                            console.log('selecting time automatically');
                            selected_time_slice = null;
                            last_time_slice = null;
                        }else{
                            selected_time_slice = last_time_slice;
                        }
                        time_state.pt = phylo_table;
                        time_state.st = sample_table;
                        time_state.tl = temp_list;
                        time_state.sid = sid;
                        time_state.mode = 'distance';
                        refresh_temporal_controls_timeline()
                        update(root);
                    }
                }
            }
        }
    }
    window.addEventListener("resize",update);
    window.addEventListener("resize",refresh_temporal_controls_timeline);
}

$(function() {
    //upper data bar for display of SID, relative abundance and species name
    $('<div id="nav_bar"></div>').appendTo('body');
    $('<span id="nav_site">phis_viz</span>').appendTo('#nav_bar');
    $('<span id="nav_data"></span>').appendTo('#nav_bar');
    $('<span id="nav_search"></span>').appendTo('#nav_bar');
    $('<span id="map_toggle"><button id="map_toggle_button">map\></button></span>').appendTo('#nav_bar');
    $('<span id="data_check"><button id="data_check_button">data\?</button></span>').appendTo('#nav_bar');

    //lower temporal bar for drilling down into a SID data set once it has already been selected
    $('<div id="temporal_bar"></div>').appendTo('body')

    //left and right graph panel containers
    $('<div id="phis_viz"></div>').appendTo('body');
    $('<span id="phis_side"></span>').appendTo('#phis_viz');
    $('<span id="map_side"></span>').appendTo('#phis_viz');
    phlyo_tree_map_graph('data/phylo_tree.csv','data/sample_data.csv','phis_side',
                         'data/site_geo.json','map_side',false);

});