

$(document).ready(function() {
        $('#refreshLogsButton').on('click', function() {
            fetchLogs();  
        });
        $('#refreshButton').on('click', function() {
            updateProgress();  
        });
    });
    




async function fetchLogs() {
    try {
        const response = await fetch('/common/fetchLogs', {
            method: 'POST', headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ task_id: taskId })});
        const data = await response.json();
        document.getElementById('log-content').textContent = data.logs;
    } catch (error) {
        console.error('Error in fetchLogs:', error);
    }
}


async function updateProgress() {
    try {
        const response = await fetch('/common/progress', {
            method: 'POST', headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ task_id: taskId })});
        const data = await response.json();
        document.getElementById('data_info').innerHTML = data[0];
    } catch (error) {
        console.error('Error in updateProgress:', error);
    }
}

    
