$(document).ready(function() {
    initializeInputSuggestions('#folderInput', 'suggestions');

    document.getElementById('selectInputDir').addEventListener('click', selectDirSteps);
    document.getElementById('executeBtn').addEventListener('click', executeSteps);
});


// selectDirSteps上传select_raw_folder，获取并创建全部脚本
// executeSteps上传task_id，执行已经构建完成的任务
// 审阅时间：20241021

async function selectDirSteps() {
    try {
        const response = await fetch('/batch_code/submit_config_script', {
            method: 'POST', headers: {'Content-Type': 'application/json'},
            body: JSON.stringify({ select_raw_folder: document.getElementById('folderInput').value})});
        const data = await response.json();
        document.getElementById('dynamicCode').textContent = data.my_code;
        hljs.highlightElement(document.getElementById('dynamicCode'));
        document.getElementById('folderInput').disabled = true;
        document.getElementById('selectInputDir').disabled = true;
        document.getElementById('executeBtn').disabled = false;
    } catch (error) {console.error('Error submitting directory:', error);}
}

async function executeSteps() {
    try {
        const response = await fetch('/common/submit_smk_job', {
            method: 'POST', headers: {'Content-Type': 'application/json'},
            body: JSON.stringify({ task_id: task_id })
        }); if (response.ok) {window.location.href = '/query/' + task_id;} 
    } catch (error) {console.error('Error submitting job:', error);}
}


