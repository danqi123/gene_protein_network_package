import os
from flask import Flask, flash, request, redirect, url_for, render_template
from plab2.startup import home_dir
from werkzeug.utils import secure_filename
from flask import session
import shutil

upload_data_path = os.path.join(str(home_dir), ".wangd0", "data", "uploads")
os.makedirs(upload_data_path, mode=0o777, exist_ok=True)

UPLOAD_FOLDER = upload_data_path
ALLOWED_EXTENSIONS = {"csv", "tsv"}

app = Flask(__name__)
app.secret_key = "someSecretKey"

app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config['MAX_CONTENT_PATH'] = 10 * 1024 * 1024  # Max 10MB


@app.template_filter()
def datetimefilter(value, format='%Y/%m/%d %H:%M'):
    """convert a datetime to a different format."""
    return value.strftime(format)


app.jinja_env.filters['datetimefilter'] = datetimefilter


@app.route("/")
def home():
    return render_template('template.html', my_string="",
                        title="Danqi's Network Analyzer",
                        hgnc_result="-",
                        ensembl_result="-",
                        uniprot_result="-",
                        hgnc_link ="",
                        uniprot_link="",
                        data_table = "")


def check_number_file():
    """Function used to get upload files and the number of file later."""
    ppi_path = []
    node_edge_path = []
    for fpath, dirname, fnames in os.walk(app.config['UPLOAD_FOLDER']):
        for file in fnames:
            if file.endswith(".csv"):
                ppi_path.append(file)
            elif file.endswith(".tsv"):
                # n_e_path = os.path.join(app.config['UPLOAD_FOLDER'], file)
                node_edge_path.append(file)
    return ppi_path, node_edge_path


@app.route("/import_option", methods=['POST'])
def import_option():
    """import uploaded files and show stats."""
    from plab2.network import Statistics
    number_ppi, number_node_edge = len(check_number_file()[0]), len(check_number_file()[1])
    if request.form["type"] == "ppi_file":
        if number_ppi != 1:
            info_string = "Please clear the contents of upload folder and upload only one PPI file."
            data_table_html = ""

        else:
            ppi_path = os.path.join(app.config['UPLOAD_FOLDER'], check_number_file()[0][0])
            session["ppi"] = ppi_path
            session["type_of_file"] = "ppi"
            info_string = "File(s) imported successfully"

            PPI_file = session.get("ppi")
            NODE_path, EDGE_path = "node_list.tsv", "edge_list.tsv"
            s = Statistics({}, None, PPI_file, None, None)
            s.write_node_list(NODE_path)
            s.write_edge_list(EDGE_path)
            s.import_graph(EDGE_path)
            data = s.summary_statistics(False)
            data.drop('Average node connectivity', inplace=True, axis=1)

            data_table_html = data.to_html(header="true", table_id="table", index=False)
            os.remove(NODE_path)
            os.remove(EDGE_path)
    elif request.form["type"] == "node_edge_file":
        if number_node_edge != 2:
            info_string = "Please clear the contents of upload folder and upload one node and one edge file."
            data_table_html = ""
        else:
            session["node_edge_1"] = os.path.join(app.config['UPLOAD_FOLDER'], check_number_file()[1][0])
            session["node_edge_2"] = os.path.join(app.config['UPLOAD_FOLDER'], check_number_file()[1][1])
            session["type_of_file"] = "node/edge"
            info_string = "File(s) imported successfully"
            try:
                NODE_path, EDGE_path = session["node_edge_1"], session["node_edge_2"]
                s = Statistics({}, None, None, NODE_path, EDGE_path)
            except:
                EDGE_path, NODE_path = session["node_edge_1"], session["node_edge_2"]
                s = Statistics({}, None, None, NODE_path, EDGE_path)
            s.import_graph(EDGE_path)
            data = s.summary_statistics(False)
            data.drop('Average node connectivity', inplace=True, axis=1)
            data_table_html = data.to_html(header="true", table_id="table", index=False)
    return render_template('template.html', my_string="",
                           title="Danqi's Network Analyzer",
                           hgnc_result="-",
                           ensembl_result="-",
                           uniprot_result="-",
                           hgnc_link="",
                           uniprot_link="",
                           reader_info= info_string,
                           data_table = data_table_html)
        # redirect("/")

@app.route("/upload")
def upload():
    return render_template('upload.html', title="Danqi's Network Analyzer")


def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS


@app.route('/uploader', methods=['GET', 'POST'])
def upload_file():
    if request.method == 'POST':
        if request.form["upload_and_clear"] == "Clear Contents":
            shutil.rmtree(app.config['UPLOAD_FOLDER'])
            os.mkdir(app.config['UPLOAD_FOLDER'])
            info = 'Uploaded files successfully removed.'
            return render_template("upload.html", title="Danqi's Network Analyzer", information = info)

        # check if the post request has the file part
        if 'file' not in request.files:
            flash('No file part')
            return redirect(request.url)
        file = request.files['file']
        # if user does not select file, browser also
        # submit an empty part without filename
        if file.filename == '':
            flash('No file selected')
            return redirect(request.url)
        if file and allowed_file(file.filename):
            filename = secure_filename(file.filename)

            if request.form["upload_and_clear"] == "Upload":
                file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
                return redirect(url_for('home'))

        return redirect(url_for('home'))

@app.route('/identifier', methods=['POST'])
def get_identifier():
    """Function for get identifier of specific symbol."""
    from plab2.startup import HGNC_data_path, UniProt_data_path, CONN_STRING
    from plab2.Utils import Profiler
    hgnc_symbol = request.form['textbox']
    p = Profiler(hgnc_symbol)
    file_path = os.path.join(HGNC_data_path, f'{hgnc_symbol}.json')
    if not os.path.isfile(file_path):
        try:
            p.request()
            hgnc_id = p.extract()["HGNC ID"]
            ensembl_id = p.extract()["Ensembl Gene ID"]
            uniprot_id = p.extract()["UniProt ID"][0]
            HGNC_root = f"http://rest.genenames.org/fetch/symbol/{hgnc_symbol}"
            uniprot_root = f"https://www.uniprot.org/uniprot/{uniprot_id}"

            return render_template('template.html', my_string=f"{hgnc_symbol}", title="Danqi's Network Analyzer",
                                   hgnc_result=hgnc_id, ensembl_result=ensembl_id,
                                   uniprot_result=uniprot_id, hgnc_link=HGNC_root, uniprot_link=uniprot_root)
        except:
            return render_template('template.html', my_string=f"No information found: {hgnc_symbol}",
                                   title="Danqi's Network Analyzer", hgnc_result="None", ensembl_result="None",
                                   uniprot_result="None", hgnc_link=f"No link found {hgnc_symbol}",
                                   uniprot_link=f"No link found {hgnc_symbol}")

    else:
        hgnc_id = p.extract()["HGNC ID"]
        ensembl_id = p.extract()["Ensembl Gene ID"]
        uniprot_id = p.extract()["UniProt ID"][0]
        HGNC_root = f"http://rest.genenames.org/fetch/symbol/{hgnc_symbol}"
        uniprot_root = f"https://www.uniprot.org/uniprot/{uniprot_id}"
        return render_template('template.html', my_string=f"{hgnc_symbol}", title="Danqi's Network Analyzer",
                               hgnc_result=hgnc_id, ensembl_result=ensembl_id,
                               uniprot_result=uniprot_id, hgnc_link=HGNC_root, uniprot_link=uniprot_root)

    # below is the version I used extract data from DATABASE:


    # # from sqlalchemy.orm import Session
    # # from sqlalchemy import create_engine
    # # from plab2.models import query_data
    # # from plab2.startup import CONN_STRING
    # # engine = create_engine(CONN_STRING, echo=True)
    # # session = Session(bind=engine)
    # # symbol = request.form['textbox']
    # #
    # if not query_data(symbol):
    #     return render_template('template.html', my_string=f"No information found {hgnc_symbol}", title="Danqi's Network Analyzer", hgnc_result="None", ensembl_result="None", uniprot_result="None", hgnc_link=f"No link found {hgnc_symbol}", uniprot_link=f"No link found {hgnc_symbol}")
    #
    # if query_data(symbol):
    #     list_info = list(query_data(symbol).values())
    #     HGNC_root = f"http://rest.genenames.org/fetch/symbol/{hgnc_symbol}"
    #     uniprot_root = f"https://www.uniprot.org/uniprot/{list_info[0][2]}"
    #     return render_template('template.html', my_string=f"{symbol}", title="Danqi's Network Analyzer", hgnc_result = "HGNC:" + str(list_info[0][0]), ensembl_result = list_info[0][1], uniprot_result = list_info[0][2], hgnc_link = HGNC_root, uniprot_link = uniprot_root)

@app.route('/short_path', methods=['POST'])
def short_path():
    """Function used to get shortest path between two HGNC symbols."""
    from plab2.network import Analyzer

    if session.get("type_of_file") == "ppi":
        PPI_file = session.get("ppi")
        node_path, edge_path = "node_list.tsv", "edge_list.tsv"
        a = Analyzer({}, None, PPI_file, None, None)
        a.write_node_list(node_path)
        a.write_edge_list(edge_path)
        a.import_graph(edge_path)
        os.remove(node_path)
        os.remove(edge_path)
        source, target = request.form["source, target"].split(",")[0], request.form["source, target"].split(",")[1]
        try:
            first_path = a.shortest_path(source, target, print_option=False)[0]
            result_sentence = f"Shortest path for {source} and {target} is:"
            return render_template('template.html', my_string="", title="Danqi's Network Analyzer", hgnc_result="-",
                                   ensembl_result="-", uniprot_result="-", hgnc_link="", uniprot_link="",
                                   result=result_sentence, reminder="", find_path=first_path)
        except:
            result_sentence = f"No paths found between {source} and {target}."
            return render_template('template.html', my_string="", title="Danqi's Network Analyzer", hgnc_result="-",
                                   ensembl_result="-", uniprot_result="-", hgnc_link="", uniprot_link="",
                                   result=result_sentence, reminder="", find_path="")

    elif session.get("type_of_file") == "node/edge":
        source, target = request.form["source, target"].split(",")[0], request.form["source, target"].split(",")[1]
        try:
            node_path, edge_path = session.get("node_edge_1"), session.get("node_edge_2")
            a = Analyzer({}, None, None, node_path, edge_path)
        except:
            edge_path, node_path = session.get("node_edge_1"), session.get("node_edge_2")
            a = Analyzer({}, None, None, node_path, edge_path)
        a.import_graph(edge_path)
        sp = a.shortest_path(source, target, print_option=False)
        if sp:
            result_sentence = f"Shortest path for {source} and {target} is:"
            return render_template('template.html', my_string="", title="Danqi's Network Analyzer", hgnc_result="-", ensembl_result="-", uniprot_result="-", hgnc_link="", uniprot_link="", result=result_sentence, reminder="", find_path=sp[0])
        else:
            result_sentence = f"No paths found between {source} and {target}."
            return render_template('template.html', my_string="", title="Danqi's Network Analyzer", hgnc_result="-", ensembl_result="-", uniprot_result="-", hgnc_link="", uniprot_link="", result=result_sentence, reminder="", find_path="")


@app.route('/plot.png', methods=['GET', 'POST'])
def plot_png():
    """Function used to plot network graph."""
    from plab2.network import Network
    import matplotlib.pyplot as plt
    import networkx as nx
    import io
    import base64
    if session.get("type_of_file") == "ppi":
        PPI_file = session.get("ppi")
        node_path, edge_path = "node_list.tsv", "edge_list.tsv"
        n = Network({}, None, PPI_file, None, None)
        n.write_node_list(node_path)
        n.write_edge_list(edge_path)
        n.import_graph(edge_path)
    elif session.get("type_of_file") == "node/edge":
        try:
            NODE_path, EDGE_path = session["node_edge_1"], session["node_edge_2"]
            n = Network({}, None, None, NODE_path, EDGE_path)
        except:
            EDGE_path, NODE_path = session["node_edge_1"], session["node_edge_2"]
            n = Network({}, None, None, NODE_path, EDGE_path)
        n.import_graph(EDGE_path)

    node_colors, edge_colors = 'red', 'black'
    plt.figure(figsize=(18, 18))
    n.graph.pos = nx.spring_layout(n.graph, k=0.06)
    nx.draw_networkx(n.graph, pos=n.graph.pos, with_labels=True,
                         node_color=node_colors,
                         edge_color=edge_colors, alpha=0.3)
    img = io.BytesIO()
    plt.savefig(img, format='png')
    plt.close()
    img.seek(0)
    plot_url = base64.b64encode(img.getvalue()).decode('utf8')
    return render_template('plot.html', plot_url=plot_url)



if __name__ == '__main__':
    flask_port = int(os.environ.get('FLASK_PORT', '5005'))
    app.run(host='0.0.0.0', port=flask_port, debug=True)
