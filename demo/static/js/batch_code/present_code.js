$(document).ready(function() {
    initializeInputSuggestions('#folderInput', 'suggestions');
});


document.addEventListener('DOMContentLoaded', function() {
    document.getElementById('selectInputDir').addEventListener('click', async function() {
        try {
            const response = await fetch('/batch_code/submit_config_script', {
                method: 'POST', headers: {'Content-Type': 'application/json'},
                body: JSON.stringify({ select_raw_folder: document.getElementById('folderInput').value, task_id: task_id})
            });
            const responseData = await response.json();
            document.getElementById('dynamicCode').textContent = responseData.my_code;
            hljs.highlightElement(document.getElementById('dynamicCode'));
            document.getElementById('folderInput').disabled = true;
            document.getElementById('selectInputDir').disabled = true;
            document.getElementById('operationTwoBtn').disabled = false;
        } catch (error) {console.error('Error submitting directory:', error);}
    });


    document.getElementById('operationTwoBtn').addEventListener('click', function() {
        window.location.href = '/query/' + task_id;
        fetch('/submit_smk_job', {
            method: 'POST', headers: {'Content-Type': 'application/json'},
            body: JSON.stringify({ task_id: task_id })
        })
    });
});

